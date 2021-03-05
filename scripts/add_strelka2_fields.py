#!/usr/bin/env python

""" A script to add standard VCF fields to Strelka2
        tumor / normal output
    Strelka2 VCF should be normalized prior to this step
    Takes input VCF and names of tumor and normal samples
    Produces modified VCF
"""

import argparse
import os
import pysam

def name_to_gt(name):
    """ Convert names like 'het' to canonical genotypes 
        'conflict' interpreted as homozygous reference
    """
    if name == 'hom':
        return (1, 1)
    elif name == 'het':
        return (0, 1)
    elif name in ('ref', 'conflict'):
        return (0, 0)
    else:
        raise ValueError('No genotype associated with %s' % name)

def get_tumor_GT(info_sgt, ref):
    """ Produce standard GT for a Strelka2 tumor sample

        Args:
            info_sgt (str): value of info['SGT'] for pysam VariantRecord object
            ref (str): value of ref for pysam VariantRecord object
            Strelka2 SGT field may be e.g. 'ref->het' or 'GG->GC'

        Returns:
            genotype (tuple): e.g. (0, 1)
    """
    sgt = info_sgt.split("->")[-1]
    try:
        tumor_gt = name_to_gt(sgt)
    except ValueError:
        if len(set(sgt)) == 2:
            tumor_gt = (0, 1)
        elif len(set(sgt)) == 1:
            if sgt[0] == ref:
                tumor_gt = (0, 0)
            else:
                tumor_gt = (1, 1)
        else:
            raise IOError('Unhandled Strelka2 SGT value %s' % sgt)
        
    return tumor_gt

def get_normal_GT(info_nt):
    """ Produce standard GT for a Strelka2 normal sample

        Args:
            info_nt: value of info['NT'] for pysam VariantRecord object

        Returns:
            genotype (tuple): e.g. (0, 1)
    """
    return name_to_gt(info_nt)

def get_AD(vcf_sample, ref, alts):
    """ Produce standard AD for a Strelka2 sample
        
        Args:
            vcf_sample (pysam VariantRecord sample object)
            ref (str): value of ref for pysam VariantRecord object
            alts (tuple): value of alts for pysam VariantRecord object
        Returns:
            allele depths (tuple): (ref, alt)
    
        From strelka2 docs:
        SNPs:
            refCounts = Sum values FORMAT column $REF + “U”
                (e.g. if REF="A" then use the value in FORMAT/AU)
            altCounts = Sum values FORMAT column $ALT + “U”
                (e.g. if ALT="T" then use the value in FORMAT/TU)
            Indels:
                tier1RefCounts = Sum values from FORMAT/TAR
                tier1AltCounts = Sume values from FORMAT/TIR
    """
    # only indels should have TIR field
    if 'TIR' in vcf_sample:
        ref_ct = sum(vcf_sample['TAR'])
        alt_ct = sum(vcf_sample['TIR'])
    else: # SNP case
        snp_ref_key = ref[0] + "U"
        snp_alt_key = alts[0][0] + "U"
        ref_ct = sum(vcf_sample[snp_ref_key])
        alt_ct = sum(vcf_sample[snp_alt_key])

    return (ref_ct, alt_ct)

def get_AF(ad):
    """ Calculate allele frequency from allele depths 
        
        Args:
             ad (tuple): (ref, alt) depths

        Returns:
             alt allele frequency (float) 
    """
    ref_counts, alt_counts = ad
    total_counts = ref_counts + alt_counts
    if not total_counts:
        return 0.0
    return alt_counts / total_counts

def create_mod_vcf(output_path, input_path, tumor_name, normal_name):
    """ Create a new VCF in which standard FORMAT fields have been
            calculated and added=based on those provided natively by Strelka2

        Args:
            output_path (str): path to the VCF being created
            input_path (str): path to the Strelka2 VCF
            tumor_name (str): name of tumor sample in input
            normal_name (str) name of normal sample in input

        Raises:
            IOError if a sample name in the input is neither tumor_name
                nor normal_name
    """

    input_vcf = pysam.VariantFile(input_path, 'r')

    input_vcf.header.formats.add('GT', '1', 'String',
            'Converted genotype added in post for compatibility')
    input_vcf.header.formats.add('AD', 'R', 'Integer',
            ('Allelic depths for the ref and alt alleles in the order listed. '
            'Added in post for compatibility'))
    input_vcf.header.formats.add('AF', 'A', 'Float',
            ('Allele frequency, as recommended by strelka2 docs: '
            '<ALT>U/<REF>U+<ALT>U (somatic snps), '
            'TIR/TIR+TAR (somatic indels)'))

    output = pysam.VariantFile(output_path, 'w', header=input_vcf.header)

    for record in input_vcf.fetch():
        for index in (0, 1):
            sample = record.samples[index]
            if sample.name == tumor_name:
                sample['GT'] = get_tumor_GT(record.info['SGT'], record.ref)
            elif sample.name == normal_name:
                sample['GT'] = get_normal_GT(record.info['NT'])
            else:
                raise IOError('VCF contains sample %s, not passed as tumor or normal'
                        % sample.name)
            sample['AD'] = get_AD(sample, record.ref, record.alts)
            sample['AF'] = get_AF(sample['AD'])

        output.write(record)

    output.close()
    input_vcf.close()

def build_output_name(inpath, output_basename):
    """Builds an VCF.GZ output filename based on the input (VCF/VCF.GZ) name
       The filename will be the input filename with '.standard' inserted before '.vcf.gz'
       
       Args:
           inpath (str): Path to the input VCF(.GZ) file
           output_basename (str): Used as first element of output filename in place of 
               first element of name of inpath filename
       Return:
           str: Filename in the format <output_basename>.<input_nameroot>.consensus.vcf.gz
    """
    basename_split = os.path.split(inpath)[-1].split('.')
    output_fields = [output_basename] + basename_split[1:-2] + ['standard'] + basename_split[-2:]
    return '.'.join(output_fields)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
            description = 'Add standard fields to Strelka2 VCF')

    parser.add_argument('--strelka2_vcf', 
            help='Normalized Strelka2 VCF')
    parser.add_argument('--tumor_name',
            help='Name of tumor sample in VCF')
    parser.add_argument('--normal_name',
            help='Name of normal sample in VCF')
    parser.add_argument('--output_basename',
            help='String to use as basename for output file [e.g.] task ID')

    args = parser.parse_args()

    # Get output VCF location
    base_dir = os.getcwd()
    output_vcf_name = build_output_name(args.strelka2_vcf, args.output_basename)
    output_vcf_path = os.path.join(base_dir, output_vcf_name)

    # Create and index the modified VCF
    create_mod_vcf(output_vcf_path, args.strelka2_vcf, args.tumor_name, args.normal_name)
    pysam.tabix_index(output_vcf_path, preset="vcf", force=True) 
