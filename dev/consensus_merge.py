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
        def strip_chr(chrom_name):
            """ Remove leading 'chr', if any, from chrom_name and return """
            if chrom_name.startswith('chr'):
                return chrom_name[3:]
            return chrom_name

        if not isinstance(other, Variant):
            return False

        if strip_chr(self.record.chrom) < strip_chr(other.record.chrom):
            return True

        if self.record.pos < other.record.pos:
            return True

        if self.record.alts[0] < other.record.alts[0]:
            return True

        return False

    def __hash__(self):
        """ Define hash for object to permit use of set() """
        return hash((self.record.chrom, self.record.pos, self.record.alts))

def write_output_header(output_vcf, contig_list):
    """ Write all header information to consensus vcf file 
        Args:
            output_vcf (pysam.VariantFile)
            contig_list (list of pysam.libcbcf.VariantContig objects)
                to be written into new header
    """
    def strip_chr(chr_name):
        """ Remove 'chr', if any, from chromosome name and return """
        if chr_name.startswith('chr'):
            return chr_name[3:]
        return chr_name

    # Only canonical chromosomes desired in consensus output
    allowed_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'M']
    
    callers = ('Strelka2', 'Mutect2', 'VarDict', 'Lancet')

    # Version is set to VCF4.2 on creation
    # Add reference? date?
    output_vcf.header.filters.add('Consensus', None, None, 'Called by two or more callers')
    output_vcf.header.filters.add('Hotspot', None, None, 'In mutational hotspot, as defined by ????')
    output_vcf.header.info.add('MQ',1,'String', 'RMS mapping quality (normal sample)') 
    output_vcf.header.formats.add('ADC', '.', 'String', 'Allele depths for %s' % ','.join(callers))
    output_vcf.header.formats.add('DPC', '.', 'String', 'Read depths for %s' % ','.join(callers))
    output_vcf.header.formats.add('AFC', '.', 'String', 'Allele frequencies for %s' % ','.join(callers))
    for contig in sorted(contig_list, key=lambda x: x.id):
        if strip_chr(contig.name) in allowed_chroms:
            output_vcf.header.contigs.add(contig.name, contig.length)

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

    return called_in

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
    write_output_header(output_vcf, strelka2_vcf.header.contigs.values())

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
    for variant in all_variants_ordered:
        hotspot = False
        seen_in = find_variant_call(variant, all_variants_dict)
        for caller_var in seen_in:
            if is_hotspot: # Need to understand from Dan's work how to check this
                hotspot = True
                output_record = build_output_record(seen_in, hotspot=True)
                output_vcf.write(output_record)
                break
        if not hotspot and len(seen_in) >= 2:
            output_record = build_output_record(seen_in)
            output_vcf.write(output_record)


