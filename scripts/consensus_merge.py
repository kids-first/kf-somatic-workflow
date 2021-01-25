#!/usr/bin/env python

""" 
From individual caller (Strelka2, Mutect2, VarDict, Lancet)
somatic vcfs, create single vcf of all variants that PASSed 
in at least 2 of the above.

Information about AD, DP, AF, GT from individual callers
should appear in the final as FORMAT fields. 

Mapping quality information is also retrieved from the CRAM
and added as an INFO field

So inputs are the four vcfs and the CRAM
All VCFs must have any 'hotspot', PON annotations added
"""
import os
import sys

import argparse
from collections import Counter
import functools
import itertools
import numpy as np
import pysam

# Ordered list of names of callers being used for consensus
CALLER_NAMES = ('Strelka2', 'Mutect2', 'VarDict', 'Lancet')
# Only canonical chromosomes desired in consensus output
# Note that these names with 'chr' prepended are also allowed
ALLOWED_CHROMS = [str(i) for i in range(0, 23)] + ['X', 'Y', 'M']
# Character for separating caller information in FORMAT fields
FORMAT_JOIN = '|'

def strip_chr(chrom_name):
    """ Remove leading 'chr', if any, from chrom_name and return """
    if chrom_name.startswith('chr'):
        return chrom_name[3:]
    return chrom_name

def stringify(iterable, sep_char=','):
    """ Turn an iterable into a delimited string with no spaces 
        delimiting character is ',' by default
    """
    return sep_char.join([str(i) for i in iterable])

class Sample(object):
    """ Object to calculate and store FORMAT information
           for samples associated with a VCF record (Variant)
    """
    @staticmethod
    def AF_from_AD(AD):
        """ Calculate allele frequency from AD FORMAT field value
        """
        ref_counts, alt_counts = AD
        total_counts = ref_counts + alt_counts
        if not total_counts:
            return 0.0
        return alt_counts / total_counts

    @staticmethod
    def name_to_gt(name):
        """ Convert names like 'het' to canonical genotypes """
        if name == 'hom':
            return '1/1'
        elif name == 'het':
            return '0/1'
        elif name in ('ref', 'conflict'):
            return '0/0'
        else:
            raise ValueError('No genotype associated with %s' % name)

    def __init__(self, name, stype, record, caller):
        """ Initialize this object and calculate 
            Args:
                name (str): the name of the sample from the VCF
                stype (str): sample may be 'tumor' or 'normal'
                record (pysam.VariantRecord): record with which this
                sample is associated
                caller (str): name of caller for this variant, 
                    must be a member of CALLER_NAMES (case insensitive)

            Raises:
                ValueError if stype is not 'tumor' or 'normal'
                IOError if name is not associated with exactly one
                    sample in record
                ValueError if caller is not in CALLER_NAMES
        """
        self.name = name
        if stype not in ('tumor', 'normal'):
            raise ValueError('Sample type must be tumor or normal, not %s' % stype)
        self.type = stype
        self.record = record

        target_sample = [self.record.samples[index] for index in range(len(self.record.samples)) 
                         if self.record.samples[index].name == self.name]
        if len(target_sample) != 1:
            raise IOError('Sample name %s not found or duplicated' % self.name)
        self.vcf_sample = target_sample[0]

        caller_test = [c.lower() for c in CALLER_NAMES]
        if caller.lower() not in caller_test:
            raise ValueError('Unknown caller %s' % caller)

        self.caller = caller.lower()

    @functools.cached_property
    def GT(self):
        """ Instructions for GT FORMAT attribute """
        if self.caller == 'strelka2':
            if self.type == 'normal':
                strelka_nt = self.record.info['NT']
                return Sample.name_to_gt(strelka_nt)
            elif self.type == 'tumor':
                # Strelka2 SGT field may be e.g. 'ref->het' or 'GG->GC'
                strelka_sgt = self.record.info['SGT'].split("->")[-1]
                try:
                    tumor_gt = Sample.name_to_gt(strelka_sgt)
                except ValueError:
                    if len(set(strelka_sgt)) == 2:
                        tumor_gt = '0/1'
                    elif len(set(strelka_sgt)) == 1:
                        if strelka_sgt[0] == self.record.ref:
                            tumor_gt = '0/0'
                        else:
                            tumor_gt = '1/1'
                    else:
                        raise IOError('Unhandled Strelka2 SGT value %s' % record.info['SGT'])
                return tumor_gt

        elif self.caller in ('mutect2', 'vardict', 'lancet'):
            return stringify(sorted(self.vcf_sample['GT']), '/')

        else:
            raise IOError('No rule to get GT for caller %s' % self.caller)

    @functools.cached_property
    def AD(self):
        """ Instructions for AD FORMAT attribute """
        if self.caller == 'strelka2':
            """ From strelka2 docs:
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
            if 'TIR' in self.vcf_sample:
                ref_ct = sum(self.vcf_sample['TAR'])
                alt_ct = sum(self.vcf_sample['TIR'])
            else: # SNP case
                snp_ref_key = self.record.ref[0] + "U"
                snp_alt_key = self.record.alts[0][0] + "U"
                ref_ct = sum(self.vcf_sample[snp_ref_key])
                alt_ct = sum(self.vcf_sample[snp_alt_key])
            
            return (ref_ct, alt_ct)

        elif self.caller in ('mutect2', 'vardict', 'lancet'):
            return self.vcf_sample['AD']

        else:
            raise IOError('No rule to get AD for caller %s' % self.caller)

    @functools.cached_property
    def AF(self):
        """ Instructions for AF FORMAT attribute """

        if self.caller == 'mutect2':
            return self.vcf_sample['AF'][0]
        elif self.caller == 'vardict':
            return self.vcf_sample['AF']
        elif self.caller == 'lancet':
            return Sample.AF_from_AD(self.vcf_sample['AD'])
        elif self.caller == 'strelka2':
            return Sample.AF_from_AD(self.AD)
        else:
            raise IOError('No rule to get AF for caller %s' % self.caller)

    @functools.cached_property
    def DP(self):
        """ Instructions for DP FORMAT attribute """
        if self.caller in ('strelka2', 'mutect2', 'vardict', 'lancet'):
            return self.vcf_sample['DP']
        else:
            raise IOError('No rule to get DP for caller %s' % self.caller)

class Variant(object):
    """ Capture relevant information and enable desired operations
            for VCF records
    """

    def __init__(self, variantrecord, caller):
        """ Initialize this object
            Args:
                variantrecord (pysam.VariantRecord)
                caller (str): the name of a caller
            Raises:
                ValueError if the record contains more than
                    one alternate allele
        """
        self.record = variantrecord
        self.caller = caller

        # VCFs should normalized prior to use with this class
        if len(self.record.alts) > 1:
            raise ValueError('Call in %s at %s:%s has multiple alternates'
                    % (caller, self.record.chrom, self.record.pos))
        
    def __eq__(self, other):
        """ Two Variants are considered equal if they occur on
               the same chromosome and at the same position and
               have the same alternate allele
        """
        if not isinstance(other, Variant):
            return False

        if (self.record.chrom == other.record.chrom
                and self.record.pos == other.record.pos
                and self.record.ref == other.record.ref
                and self.record.alts == other.record.alts):
            return True

        return False

    def __lt__(self, other):
        """ Define whether this object is less than another Variant
            First based on the chromosome name, stripped of leading 'chr'
            Then based on chromosomal coordinate, as 'int'
            Then alphabetically by first alternate allele
        """
        def order_chroms(chrom1, chrom2):
            """ Establish correct sort order for canonical human chromosomes
                Returns True if chrom1 < chrom2, else False
            """
            chr1 = strip_chr(chrom1)
            chr2 = strip_chr(chrom2)
            
            for name in chr1, chr2:
                if name not in ALLOWED_CHROMS:
                    raise ValueError('Uncanonical chromosome %s uncomparable' % name)

            letter_mappings = {'X': 23, 'Y': 24, 'M': 25}
            try:
                chr1_num = letter_mappings[chr1]
            except KeyError:
                chr1_num = int(chr1)
            try:
                chr2_num = letter_mappings[chr2]
            except KeyError:
                chr2_num = int(chr2)

            return (chr1_num < chr2_num)


        if not isinstance(other, Variant):
            raise ValueError('Cannot compare Variant %s to non-Variant %s'
                             % (self, other))

        if order_chroms(self.record.chrom, other.record.chrom):
            return True
        elif self.record.chrom != other.record.chrom:
            return False
        else:
            if self.record.pos < other.record.pos:
                return True
            elif self.record.pos != other.record.pos:
                return False
            else:
                if self.record.ref < other.record.ref:
                    return True
                elif self.record.ref != other.record.ref:
                    return False
                else:
                    if self.record.alts[0] < other.record.alts[0]:
                        return True 
                    else:
                        return False

        return False

    def __hash__(self):
        """ Define hash for object to permit use of set() """
        return hash((self.record.chrom, self.record.pos, self.record.ref, self.record.alts))

    def __repr__(self):
        """ How the object looks when printed """
        return ' '.join([self.record.chrom, str(self.record.start), 
                         self.record.ref, str(self.record.alts), self.caller])

def write_output_header(output_vcf, sample_list, contig_list):
    """ Write all header information to consensus vcf file 
        Args:
            output_vcf (pysam.VariantFile)
            sample_list (pysam.libcbcf.VariantHeaderSamples): samples to write into header
            contig_list (list of pysam.libcbcf.VariantContig objects):
                names of contigs to write into header
    """
    # Version is set to VCF4.2 on creation
    # Add reference? date?
    output_vcf.header.info.add('MQ',1,'Integer', 
            'RMS mapping quality (normal sample)')
    output_vcf.header.info.add('MQ0',1,'Integer', 
            'Number of MAPQ=0 reads (normal sample)')
    output_vcf.header.info.add('CAL','.','String',
            'List of callers making this call')
    output_vcf.header.info.add('Hotspot', 0, 'Flag',
            'In mutational hotspot, as defined by ????')
    output_vcf.header.formats.add('GT', '1', 'String', 
            'Consensus genotype')
    output_vcf.header.formats.add('AD', 'R', 'Integer',
            'Consensus depths for the ref and alt alleles in the order listed')
    output_vcf.header.formats.add('AF', 'A', 'Float',
            'Consensus allele frequency')
    output_vcf.header.formats.add('DP', '1', 'Integer',
            'Consensus depth')
    output_vcf.header.formats.add('GTC', '.', 'String', 
            'Genotypes for %s' % FORMAT_JOIN.join(CALLER_NAMES))
    output_vcf.header.formats.add('GT_STATUS', '.', 'String', 
            ("Degree of unanimity of genotype: 'unanimous', 'majority', or 'deadlock' if" 
            " equally supported by individual callers"))
    output_vcf.header.formats.add('ADC', '.', 'String', 
            'Allele depths for %s' % FORMAT_JOIN.join(CALLER_NAMES))
    output_vcf.header.formats.add('DPC', '.', 'String', 
            'Read depths for %s' % FORMAT_JOIN.join(CALLER_NAMES))
    output_vcf.header.formats.add('AFC', '.', 'String', 
            'Allele frequencies for %s' % FORMAT_JOIN.join(CALLER_NAMES))
    for contig in sorted(contig_list, key=lambda x: x.id):
        if strip_chr(contig.name) in ALLOWED_CHROMS:
            output_vcf.header.contigs.add(contig.name, contig.length)
    for sample in sample_list:
        output_vcf.header.add_sample(sample)

def get_all_variants(caller_name, pysam_vcf):
    """ Ingest all records from a pysam VCF object

        Args:
            caller_name (str)
            pysam_vcf (pysam.VariantFile)

        Returns:
            dict: with key caller_name and value list of Variant objects.
                'caller' attribute of Variants will be caller_name.
    """

    all_records = sorted([Variant(rec, caller_name) for rec in pysam_vcf.fetch()])
    return {caller_name: all_records}

def find_variant_callers(variant_list, all_variants_dict):
    """ Find callers that have called each of a list of Variants
        Args:
            variant_list (list of Variants): variants being searched for
            all_variants_dict (dict): associates caller names
                with lists of Variants
        Return:
            list of lists of single-caller Variants

        NOTE: This requires all lists to be sorted but does not check!
    """
    single_caller_variants = []
    for var in variant_list:
        sublist = []
        for caller, varlist in all_variants_dict.items():
           if not varlist:
               continue
           if varlist[0] == var:
               sublist.append(varlist.pop(0))
        single_caller_variants.append(sublist)

    return single_caller_variants

def get_gt_consensus(gt_list):
    """ Determine level of genotype consensus from single-caller genotypes 
        Args:
            gt_list (list): Strings in '0/0' format representing 
                single-caller genotypes

        Returns: tuple of (consensus genotype, consensus tag)
            consensus genotype will be '0/1' in event of tie
            consensus tag may be 'unanimous', 'majority', 'deadlock'
    """
    conflict_gt_tag = '0/1'
    unan_tag, majo_tag, tie_tag = 'unanimous', 'majority', 'deadlock'

    counts = Counter(gt_list)

    if '.' in counts:
        counts.pop('.')

    if len(counts) == 1:
        return list(counts)[0], unan_tag

    largest_count = max(counts.values())
    has_max_count = {k: v for k, v in counts.items() if v == largest_count}
    if len(has_max_count) == 1:
        return list(has_max_count)[0], majo_tag
    else:
        return conflict_gt_tag, tie_tag

def get_ad_consensus(ad_list):
    """ Get mean value of AD across callers 
        Args:
           ad_list (list): list of AD (comma-separated strings)
        Return:
           (tuple): mean depths
    """
    ref_depth = int(np.mean([int(f.split(',')[0]) for f in ad_list if f != '.']))
    alt_depth = int(np.mean([int(f.split(',')[1]) for f in ad_list if f != '.']))
    return (ref_depth, alt_depth)

def get_af_consensus(af_list):
    """ Get mean value of AF across callers
        Args:
            af_list (list): list of AF (strings)
        Return:
            (float): mean AF
    """
    return np.mean([float(f) for f in af_list if f != '.'])

def get_dp_consensus(dp_list):
    """ Get mean value of DP across callers
        Args:
            dp_list (list): list of DP (strings)
        Return:
            (int): mean DP
    """
    return int(np.mean([int(f) for f in dp_list if f != '.']))

def get_mapq(cram, chrom, pos):
    """ Get RMS mapping quality and number of MAPQ=0 reads at locus

        Args:
            cram (pysam.AlignmentFile): open in 'r' mode
            chrom (str): locus chromosome name
            pos (int): 1-based locus coordinate

        Returns:
            tuple of ints RMS mapq and # of mapq reads at locus
    """
    mapq = []
    mq0 = 0
#    print(chrom, pos-1, pos)

    aligned_reads = cram.fetch(chrom, pos-1, pos)

    for read in aligned_reads:
        mapq.append(read.mapq)
        if read.mapq == 0:
            mq0 += 1

    if not mapq:
        return 0, 0

    rms_mapq = np.sqrt(np.mean(np.array(mapq)**2))
    return int(rms_mapq), mq0

def build_output_record(single_caller_variants, output_vcf, normal_cram, hotspot=False):
    """ Create a single VCF record for the consensus output
            by interrogating the single-caller records for that variant
            and add that record the consensus VCF

        Args: 
            single_caller_variants (list of Variants): list of Variants
                representing a single call in one or more callers
            output_vcf (pysam.VariantFile): consensus file, open in 'w' mode
            normal_cram (pysam.Alignmentfile): Normal CRAM for this sample,
                open in 'r' mode
            hotspot (bool): whether this variant is in a mutational hotspot
                (default False)
    """
    variant_lookup = {v.caller: v for v in single_caller_variants}

    output_record = output_vcf.new_record()
    # Set consistent attributes
    liftover_record = single_caller_variants[0].record
    output_record.chrom = liftover_record.chrom
    output_record.ref = liftover_record.ref
    output_record.alts = liftover_record.alts
    output_record.id = liftover_record.id
    output_record.start = liftover_record.start
    output_record.stop = liftover_record.stop

    # For each caller, get information required for format fields for this variant
    joint_tags = ('GT', 'AD', 'AF', 'DP')
    no_call_info = {tag: '.' for tag in joint_tags} 
    allele_information = {caller.lower(): {} for caller in CALLER_NAMES}

    for variant in single_caller_variants:
        normal_sample, tumor_sample = variant.record.samples
        variant.normal = Sample(normal_sample, 'normal', variant.record, variant.caller)
        variant.tumor = Sample(tumor_sample, 'tumor', variant.record, variant.caller)

    for index in (0, 1):
        normal = 0
        tumor = 1
        GT_list = []
        AD_list = []
        AF_list = []
        DP_list = []
        for caller in CALLER_NAMES:
            try:
                targ_var = variant_lookup[caller.lower()]
            except KeyError:
                for flist in (GT_list, AD_list, AF_list, DP_list):
                    flist.append('.')
                continue
            if index == normal:
                targ_sample = targ_var.normal
            elif index == tumor:
                targ_sample = targ_var.tumor
                
            GT_list.append(str(targ_sample.GT))
            AD_list.append(stringify(targ_sample.AD))
            AF_list.append('{:0.4f}'.format(targ_sample.AF))
            DP_list.append(str(targ_sample.DP))

        GT = stringify(GT_list, FORMAT_JOIN)
        AD = stringify(AD_list, FORMAT_JOIN)
        AF = stringify(AF_list, FORMAT_JOIN)
        DP = stringify(DP_list, FORMAT_JOIN)

        output_record.samples[index]['GTC'] = GT
        output_record.samples[index]['ADC'] = AD
        output_record.samples[index]['AFC'] = AF
        output_record.samples[index]['DPC'] = DP

        consensus_gt, gt_tag = get_gt_consensus(GT_list)

#        print(consensus_gt, tuple(consensus_gt.split('/')))
        output_record.samples[index]['GT'] = tuple([int(i) for i in consensus_gt.split('/')])
        output_record.samples[index]['GT_STATUS'] = gt_tag
        output_record.samples[index]['AD'] = get_ad_consensus(AD_list)
        output_record.samples[index]['AF'] = get_af_consensus(AF_list)
        output_record.samples[index]['DP'] = get_dp_consensus(DP_list)

    chrom = single_caller_variants[0].record.chrom
    pos = single_caller_variants[0].record.pos
    mapq, mq0 = get_mapq(normal_cram, chrom, pos)

    output_record.info['MQ'] = mapq
    output_record.info['MQ0'] = mq0
    output_record.info['CAL'] = ','.join([c for c in CALLER_NAMES if c.lower() in variant_lookup])

    if hotspot:
        output_record.info['Hotspot'] = True

    output_record.filter.add('PASS')

    output_vcf.write(output_record)

def build_output_name(inpath, outbase):
    """Builds an VCF.GZ output filename based on the input (VCF/VCF.GZ) name and given output basename
       Args:
           inpath (str): Path to the input VCF(.GZ) file
           outbase (str): Given output basename string
       Return:
           str: Filename in the format <output_basename>.<input_nameroot>.consensus.vcf.gz
    """
    basename_split = os.path.split(inpath)[-1].split('.')
    return '.'.join([outbase, basename_split[0] + '.consensus.vcf.gz'])

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--strelka2_vcf', 
                        help='Strelka2 VCF with MNPs fixed')
    parser.add_argument('--mutect2_vcf',
                        help='Mutect2 VCF')
    parser.add_argument('--lancet_vcf', 
                        help='Lancet VCF')
    parser.add_argument('--vardict_vcf',
                        help='VarDict VCF')
    parser.add_argument('--cram', help='CRAM file')
    parser.add_argument('--output_basename',
                        help='String to use as basename for output file')

    args = parser.parse_args()

    # Create VCF and CRAM pysam objects
    strelka2_vcf = pysam.VariantFile(args.strelka2_vcf, 'r')
    mutect2_vcf = pysam.VariantFile(args.mutect2_vcf, 'r')
    lancet_vcf = pysam.VariantFile(args.lancet_vcf, 'r')
    vardict_vcf = pysam.VariantFile(args.vardict_vcf, 'r')

    normal_cram = pysam.AlignmentFile(args.cram, 'rc') #, reference_filename="/home/ubuntu/volume/ref/Homo_sapiens_assembly38.fasta.fai")
    
    # Create output vcf
    base_dir = os.path.split(os.path.abspath(args.strelka2_vcf))[0]
    output_vcf_name = build_output_name(args.vardict_vcf, args.output_basename)
    output_vcf_path = os.path.join(base_dir, output_vcf_name)
    # print(output_vcf_path)
    output_vcf = pysam.VariantFile(output_vcf_path, 'w')
    write_output_header(output_vcf, strelka2_vcf.header.samples, 
                        strelka2_vcf.header.contigs.values())

    all_variants_dict = {}

    # Read in all variant records
    for caller, vcf in [('strelka2', strelka2_vcf),
                        ('mutect2', mutect2_vcf),
                        ('lancet', lancet_vcf),
                        ('vardict', vardict_vcf)]:
        all_variants_dict.update(get_all_variants(caller, vcf))

    # Get single ordered list of all variants
    all_variants_list = list(itertools.chain.from_iterable(all_variants_dict.values()))
    all_variants_ordered = sorted(list(set(all_variants_list)))
    single_caller_variants = find_variant_callers(all_variants_ordered, all_variants_dict)
 
    # Identify variants meeting consensus criteria
    # Current criteria are being in a mutational hotspot 
    #    or having been called by 2 or more callers
    for index, varlist in enumerate(single_caller_variants):
        if not varlist:
            break

        hotspot = False

        for caller_var in varlist:
            try:
                if caller_var.record.info['HotSpotAllele']:
                    hotspot = True
                    build_output_record(varlist, output_vcf, hotspot=True)
                    break
            except KeyError:
                continue
 
        if not hotspot and len(varlist) >= 2:
            build_output_record(varlist, output_vcf, normal_cram)

    normal_cram.close()
    output_vcf.close()
