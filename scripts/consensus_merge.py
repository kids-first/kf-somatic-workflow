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
from collections import Counter
import functools
import itertools
import os

import multiprocessing as mp
import numpy as np
import pysam

# Ordered list of names of callers being used for consensus
CALLER_NAMES = ('Strelka2', 'Mutect2', 'VarDict', 'Lancet')
# Character for separating caller information in FORMAT fields
FORMAT_JOIN = '|'

def get_range(iterable):
    """ Get the difference between the minimum and maximum values of
           the float-castable items in iterable
        Returns 0 if the number of such values is less than 2
    """
    floats = []
    for item in iterable:
        try:
            floats.append(float(item))
        except ValueError:
            continue
    if len(floats) < 2:
        return 0.0
    sorted_floats = sorted(floats)
    return sorted_floats[-1] - sorted_floats[0]

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

    def __init__(self, name, record, caller):
        """ Initialize this object and calculate
            Args:
                name (str): the name of the sample from the VCF
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
        self.record = record

        target_sample = [self.record.samples[index] for index in range(len(self.record.samples))
                         if self.record.samples[index].name == self.name]
        if len(target_sample) != 1:
            raise IOError('Sample name %s not found or duplicated in %s' % (self.name, caller))
        self.vcf_sample = target_sample[0]

        if caller.lower() not in [c.lower() for c in CALLER_NAMES]:
            raise ValueError('Unknown caller %s' % caller)

        self.caller = caller.lower()

    @functools.cached_property
    def GT(self):
        """ Instructions for GT FORMAT attribute """

        if self.caller in ('mutect2', 'vardict', 'lancet', 'strelka2'):
            return stringify(sorted(self.vcf_sample['GT']), '/')

        raise IOError('No rule to get GT for caller %s' % self.caller)

    @functools.cached_property
    def AD(self):
        """ Instructions for AD FORMAT attribute """

        if self.caller in ('mutect2', 'vardict', 'lancet', 'strelka2'):
            return self.vcf_sample['AD']

        raise IOError('No rule to get AD for caller %s' % self.caller)

    @functools.cached_property
    def AF(self):
        """ Instructions for AF FORMAT attribute """

        if self.caller in ('mutect2', 'strelka2'):
            return self.vcf_sample['AF'][0]
        if self.caller == 'vardict':
            return self.vcf_sample['AF']
        if self.caller == 'lancet':
            return Sample.AF_from_AD(self.vcf_sample['AD'])

        raise IOError('No rule to get AF for caller %s' % self.caller)

    @functools.cached_property
    def DP(self):
        """ Instructions for DP FORMAT attribute """
        if self.caller in ('strelka2', 'mutect2', 'vardict', 'lancet'):
            return self.vcf_sample['DP']

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
        def order_chroms(chr1, chr2):
            """ Establish correct sort order for canonical human chromosomes
                Returns True if chrom1 < chrom2, else False
            """
            return ALLOWED_CHROMS.index(chr1) < ALLOWED_CHROMS.index(chr2)

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

def write_output_header(output_vcf, sample_list, contig_list, hotspot_source=None):
    """ Write all header information to consensus vcf file
        Args:
            output_vcf (pysam.VariantFile)
            sample_list (pysam.libcbcf.VariantHeaderSamples): samples to write into header
            contig_list (list of pysam.libcbcf.VariantContig objects):
                names of contigs to write into header
            hotspot_source (str): information about source of 'HotSpotAllele' designated regions
                for inclusion in header; default None
    """
    # Version is set to VCF4.2 on creation
    # Add reference? date?

    # Hotspot source as of 01/2021 was
    # Memorial Sloan Kettering Cancer Center
    # based on Chang et al. 2017; see https://www.cancerhotspots.org
    if hotspot_source:
        hotspot_string = ' as defined by %s' % hotspot_source
    else:
        hotspot_string = ''

    output_vcf.header.info.add('MQ',1,'Integer',
            'RMS mapping quality (normal sample)')
    output_vcf.header.info.add('MQ0',1,'Integer',
            'Number of MAPQ=0 reads (normal sample)')
    output_vcf.header.info.add('CAL','.','String',
            'List of callers making this call')
    output_vcf.header.info.add("HotSpotAllele",'A','Integer',
            'Included by exception to consensus rule due to hotspot status%s' % hotspot_string)
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
    output_vcf.header.formats.add('ADR', 'R', 'Integer',
            'Difference between highest and lowest AD')
    output_vcf.header.formats.add('AFR', 'A', 'Float',
            'Difference between highest and lowest AF')
    output_vcf.header.formats.add('DPR', '1', 'Integer',
            'Difference between highest and lowest DP')
    for contig in sorted(contig_list, key=lambda x: x.id):
        if contig.name in ALLOWED_CHROMS:
            output_vcf.header.contigs.add(contig.name, contig.length)
    for sample in sample_list:
        output_vcf.header.add_sample(sample)

def get_all_variants(caller_name, pysam_vcf):
    """ Ingest all records from a pysam VCF object
        Excludes records associated with chromosomes not allowed by user

        Args:
            caller_name (str)
            pysam_vcf (pysam.VariantFile)

        Returns:
            dict: with key caller_name and value list of Variant objects.
                'caller' attribute of Variants will be caller_name.
    """

    all_records = sorted([Variant(rec, caller_name) for rec in pysam_vcf.fetch() 
                          if rec.chrom in ALLOWED_CHROMS])
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
           # remove other copies of the same variant
           if varlist:
               while varlist[0] == var:
                   varlist.pop(0)
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

def get_mapq(cram_path, chrom, pos, reference=None):
    """ Get RMS mapping quality and number of MAPQ=0 reads at locus

        Args:
            cram_path (str): path to CRAM or BAM file
            chrom (str): locus chromosome name
            pos (int): 1-based locus coordinate
            reference: path to reference FASTA, required for CRAM
                (default None)

        Raises:
            IOError if CRAM provided without reference
                    if cram_path has neither cram nor bam extension

        Return:
            tuple of ints RMS mapq and # of mapq reads at locus
    """
    if cram_path.endswith('cram'):
        if not reference:
            raise IOError('Must provide reference with CRAM')
        cram = pysam.AlignmentFile(cram_path, 'rc', reference_filename=reference)
    elif cram_path.endswith('bam'):
        cram = pysam.AlignmentFile(cram_path, 'rb')
    else:
        raise IOError('File provided to "cram" argument must have .cram or .bam extension')

    mapq = []
    mq0 = 0

    aligned_reads = cram.fetch(chrom, pos-1, pos)

    for read in aligned_reads:
        mapq.append(read.mapq)
        if read.mapq == 0:
            mq0 += 1

    if not mapq:
        return 0, 0

    rms_mapq = np.sqrt(np.mean(np.array(mapq)**2))
    return int(rms_mapq), mq0

def build_output_record(single_caller_variants, output_vcf, sample_names, hotspot=False):
    """ Create a single VCF record for the consensus output
            by interrogating the single-caller records for that variant
            and add that record the consensus VCF

        Args:
            single_caller_variants (list of Variants): list of Variants
                representing a single call in one or more callers
            output_vcf (pysam.VariantFile): consensus file, open in 'w' mode
            sample_names (list of str): names of samples from a single input vcf
            hotspot (bool): whether this variant is in a mutational hotspot
                (default False)
    """
    output_record = output_vcf.new_record()
    # Set consistent attributes
    liftover_record = single_caller_variants[0].record
    output_record.chrom = liftover_record.chrom
    output_record.alleles = liftover_record.alleles
    output_record.id = liftover_record.id
    output_record.start = liftover_record.start
    output_record.stop = liftover_record.stop

    # For each caller, get information required for format fields for this variant
    for variant in single_caller_variants:
        variant.samples = [Sample(name, variant.record, variant.caller)
                for name in sample_names]

    variant_lookup = {v.caller: v for v in single_caller_variants}

    for index in range(len(sample_names)):
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
            targ_sample = targ_var.samples[index]

            GT_list.append(str(targ_sample.GT))
            AD_list.append(stringify(targ_sample.AD))
            AF_list.append('{:0.4f}'.format(targ_sample.AF))
            DP_list.append(targ_sample.DP)

        GT = stringify(GT_list, FORMAT_JOIN)
        AD = stringify(AD_list, FORMAT_JOIN)
        AF = stringify(AF_list, FORMAT_JOIN)
        DP = stringify(DP_list, FORMAT_JOIN)

        output_record.samples[index]['GTC'] = GT
        output_record.samples[index]['ADC'] = AD
        output_record.samples[index]['AFC'] = AF
        output_record.samples[index]['DPC'] = DP

        consensus_gt, gt_tag = get_gt_consensus(GT_list)
        dp_range = get_range(DP_list)
        af_range = get_range(AF_list)
        ad_range_1 = get_range([item.split(',')[0] for item in AD_list])
        ad_range_2 = get_range([item.split(',')[-1] for item in AD_list])

        output_record.samples[index]['GT'] = tuple([int(i) for i in consensus_gt.split('/')])
        output_record.samples[index]['GT_STATUS'] = gt_tag
        output_record.samples[index]['AD'] = get_ad_consensus(AD_list)
        output_record.samples[index]['AF'] = get_af_consensus(AF_list)
        output_record.samples[index]['DP'] = get_dp_consensus(DP_list)
        output_record.samples[index]['ADR'] = (ad_range_1, ad_range_2)
        output_record.samples[index]['AFR'] = af_range
        output_record.samples[index]['DPR'] = dp_range

    output_record.info['CAL'] = ','.join([c for c in CALLER_NAMES if c.lower() in variant_lookup])

    output_record.info['HotSpotAllele'] = (int(hotspot),)

    output_record.filter.add('PASS')

    return output_record

def allowed_contigs_from_bed(bed_path, known_contigs):
    """ Gets list of contigs from BED file
        File must be tab-delimited with contig names in first field
        Order of contigs will be same as in list of known contigs

        Args:
            bed_path (str): path to BED file
            known_contigs (list): valid contig names (strings), presumably from VCF header

        Raises:
            IOError if bed_path does not point to a file
                    if no contig names are identified
                    if any contig name is not in known_contigs

        Return:
            list: list of contig names (strings)
    """

    if not os.path.isfile(bed_path):
        raise IOError('Contig BED file not found at %s' % bed_path)

    contig_names = []
    with open(bed_path, 'r') as bed:
        for line in bed:
            fields = line.rstrip().split('\t')
            name = fields[0]
            if name in contig_names:
                continue
            if name not in known_contigs:
                raise IOError('Contig %s not known from e.g. VCF header' % name)
            contig_names.append(name)

    if not contig_names:
        raise IOError('No contig names identified in %s' % bed_path)

    return sorted(contig_names, key=lambda name: known_contigs.index(name))

def build_output_name(inpath, output_basename):
    """Builds an VCF.GZ output filename based on the input (VCF/VCF.GZ) name
       The filename will be the input filename with '.consensus' inserted before '.vcf.gz'
           and a change of basename to the provided value

       Args:
           inpath (str): Path to the input VCF(.GZ) file
           output_basename (str): Used as first element of output filename in place of 
               first element of name of inpath filename
       Return:
           str: Filename in the format <output_basename>.<input_nameroot>.consensus.vcf.gz
    """
    basename_split = os.path.split(inpath)[-1].split('.')
    output_fields = [output_basename] + basename_split[1:-2] + ['consensus'] + basename_split[-2:]
    return '.'.join(output_fields)

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
    parser.add_argument('--cram', help='CRAM or BAM file')
    parser.add_argument('--reference',
                        help='Path to FASTA to which CRAM is aligned')
    parser.add_argument('--ncallers', type=int, default=2,
                        help='Number of callers required for consensus [2]')
    parser.add_argument('--output_basename',
                        help='String to use as basename for output file')
    parser.add_argument('--hotspot_source',
                        help='Optional source of "hotspot" designated regions')
    parser.add_argument('--contig_bed',
                        help='Optional source of names of contigs to be included in output')

    args = parser.parse_args()

    if args.ncallers < 1:
        raise ValueError('ncallers must be positive')
    if args.ncallers > len(CALLER_NAMES):
        raise ValueError('ncallers cannot be greater than number of known callers (%s)' 
                % len(CALLER_NAMES))

    # Create VCF and CRAM pysam objects
    strelka2_vcf = pysam.VariantFile(args.strelka2_vcf, 'r')
    mutect2_vcf = pysam.VariantFile(args.mutect2_vcf, 'r')
    lancet_vcf = pysam.VariantFile(args.lancet_vcf, 'r')
    vardict_vcf = pysam.VariantFile(args.vardict_vcf, 'r')

    # Get ordered names of contigs desired in output VCF and store in global variable
    known_contigs = list(strelka2_vcf.header.contigs)
    if args.contig_bed:
        allowed_contigs = allowed_contigs_from_bed(args.contig_bed, known_contigs)
    else:
        allowed_contigs = known_contigs
    global ALLOWED_CHROMS
    ALLOWED_CHROMS = allowed_contigs

    # Create output vcf
    base_dir = os.getcwd()
    output_vcf_name = build_output_name(args.vardict_vcf, args.output_basename)
    output_vcf_path = os.path.join(base_dir, output_vcf_name)
    output_vcf = pysam.VariantFile(output_vcf_path, 'w')
    write_output_header(output_vcf, strelka2_vcf.header.samples,
                        strelka2_vcf.header.contigs.values(), args.hotspot_source)

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
    sample_names = list(strelka2_vcf.header.samples)
    all_output_records = []
    for index, varlist in enumerate(single_caller_variants):
        if not varlist:
            if not index:
                raise IOError('First variant record empty')
            raise IOError('Empty record following %s' % varlist[index-1][0])

        hotspot = False
        record = None

        for caller_var in varlist:
            try:
                # HotSpotAllele INFO field is stored as length-1 tuple 
                if caller_var.record.info['HotSpotAllele'][0]:
                    hotspot = True
                    record = build_output_record(varlist, output_vcf, sample_names,
                            hotspot=True)
                    break
            except KeyError:
                continue

        if not hotspot and len(varlist) >= args.ncallers:
            record = build_output_record(varlist, output_vcf, sample_names)

        if record:
            all_output_records.append(record)

    pool = mp.Pool()
    mapq_info = pool.starmap(get_mapq, ((args.cram, rec.chrom, rec.pos, args.reference)
            for rec in all_output_records))
    pool.close()
    pool.join()

    for record, (mapq, mq0) in zip(all_output_records, mapq_info):
        record.info['MQ'] = mapq
        record.info['MQ0'] = mq0
        output_vcf.write(record)

    output_vcf.close()

    pysam.tabix_index(output_vcf_path, preset="vcf", force=True)
