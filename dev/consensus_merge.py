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
import argparse
import itertools
import os
import pysam
import sys

# Ordered list of names of callers being used for consensus
CALLER_NAMES = ('Strelka2', 'Mutect2', 'VarDict', 'Lancet')
# Only canonical chromosomes desired in consensus output
# Note that these names with 'chr' prepended are also allowed
ALLOWED_CHROMS = [str(i) for i in range(0, 23)] + ['X', 'Y', 'M']

def strip_chr(chrom_name):
    """ Remove leading 'chr', if any, from chrom_name and return """
    if chrom_name.startswith('chr'):
        return chrom_name[3:]
    return chrom_name

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
                if self.record.alts[0] < other.record.alts[0]:
                    return True
                else:
                    return False

        return False

    def __hash__(self):
        """ Define hash for object to permit use of set() """
        return hash((self.record.chrom, self.record.pos, self.record.alts))

    def __repr__(self):
        """ How the object looks when printed """
        return ' '.join([self.record.chrom, str(self.record.start), 
                         str(self.record.alts), self.caller])

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
    output_vcf.header.filters.add('Consensus', None, None, 'Called by two or more callers')
    output_vcf.header.filters.add('Hotspot', None, None, 'In mutational hotspot, as defined by ????')
    output_vcf.header.info.add('MQ',1,'String', 'RMS mapping quality (normal sample)') 
    output_vcf.header.formats.add('ADC', '.', 'String', 'Allele depths for %s' % ','.join(CALLER_NAMES))
    output_vcf.header.formats.add('DPC', '.', 'String', 'Read depths for %s' % ','.join(CALLER_NAMES))
    output_vcf.header.formats.add('AFC', '.', 'String', 'Allele frequencies for %s' % ','.join(CALLER_NAMES))
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

    all_records = [Variant(rec, caller_name) for rec in pysam_vcf.fetch()]
    return {caller_name: all_records}

def find_variant_call(variant, all_variants_dict):
    """ Find callers that have called a particular variant
        Args:
            variant (Variant): variant being searched for
            all_variants_dict (dict): associates caller names
                with lists of Variants
        Return:
            list: All Variant objects matching the searched-for variant
    """
    called_in = []
    for caller_name, var_list in all_variants_dict.items():
        for var in var_list:
            if variant == var:
                called_in.append(var)
                break

    return called_in

def build_output_record(single_caller_variants, output_vcf, hotspot=False):
    """ Create a single VCF record for the consensus output
            by interrogating the single-caller records for that variant
            and add that record the consensus VCF

        Args: 
            single_caller_variants (list of Variants): list of Variants
                representing a single call in one or more callers
            output_vcf (pysam.VariantFile): consensus file, open in 'w' mode
            hotspot (bool): whether this variant is in a mutational hotspot
                (default False)
    """
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
    joint_tags = ('AD', 'AF', 'DP')
    no_call_info = {tag: '.' for tag in joint_tags} 
    allele_information = {}
    for caller_name in CALLER_NAMES:
        for variant in single_caller_variants:
            if variant.caller == caller_name.lower():
                # allele_information[variant.caller] = get_allele_info(variant)
                allele_information[variant.caller] = no_call_info
                break
        else:
            allele_information[caller_name.lower()] = no_call_info

    output_record.info['MQ'] = '10'

    """
    for tag in joint_tags:
        format_field = '/'.join([allele_information[caller.lower()][tag] 
            for caller in CALLER_NAMES])
        output_record.format[tag] = format_field
    """
    output_vcf.write(output_record)

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

#    cram = pysam.AlignmentFile(args.cram, 'rc')
    
    # Create output vcf
    base_dir = os.path.split(os.path.abspath(args.strelka2_vcf))[0]
    output_vcf_name = args.output_basename + '.vcf.gz'
    output_vcf_path = os.path.join(base_dir, output_vcf_name)
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

    # Identify variants meeting consensus criteria
    # Current criteria are being in a mutational hotspot 
    #    or having been called by 2 or more callers
    for variant in all_variants_ordered[:100]:
        hotspot = False
        seen_in = find_variant_call(variant, all_variants_dict)
        for caller_var in seen_in:
            pass
            """
            if caller_var.record.info['HotSpotAllele']:
                hotspot = True
                build_output_record(seen_in, output_vcf, hotspot=True)
                break
            """
        if not hotspot and len(seen_in) >= 2:
            build_output_record(seen_in, output_vcf)

    output_vcf.close()
