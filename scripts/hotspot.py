#!/usr/bin/env python3

import argparse
import csv
import os
import pysam
import gzip
import re

def parse_command():
    """Function to parse the command line input
    Args:
        None
    Return:
        Namespace: returns the args as a standard Namespace object
    """
    parser = argparse.ArgumentParser(description='Flag HotSpots in VCF using one or more Hotspot File Inputs')
    required_options = parser.add_argument_group('required arguments')
    hotspot_options = parser.add_argument_group('hotspot file input arguments')
    advanced_options = parser.add_argument_group('advanced user arguments')
    hotspot_options.add_argument('--genomic_hotspots',
        nargs='*',
        help="Tab-delimited BED-formatted (chrom, zero-based chromStart, non-inclusive chromEnd) file(s) containing hg38 genomic positions corresponding to hotspots")
    hotspot_options.add_argument('--protein_snvs',
        nargs='*',
        help="Column name-labeled, tab-delimited file(s) containing HUGO-formatted protein names and VEP-formatted positions <start_aa_pos> for SNV hotspots")
    hotspot_options.add_argument('--protein_indels',
        nargs='*',
        help="Column name-labeled, tab-delimited file(s) containing HUGO-formatted protein names and VEP-formatted positions <start_aa_pos>-<end_aa_pos> for INDEL hotspots")
    required_options.add_argument('--vcf',
        required=True,
        help='Input VCF or VCF.GZ file')
    advanced_options.add_argument('--protein_colname',
        default='Hugo_Symbol',
        metavar='Hugo_Symbol',
        help='Overrides the column name in the protein_hotspots file(s) that contains the protein name information')
    advanced_options.add_argument('--position_colname',
        default='Amino_Acid_Position',
        metavar='Amino_Acid_Position',
        help='Overrides the column name in the protein_hotspots file(s) that contains the HGVSp short information')
    advanced_options.add_argument('--csq_field',
        default='SYMBOL',
        metavar='SYMBOL',
        help='Overrides the name of the CSQ field that matches the values found in the protein_colname of the protein_hotspots input(s)')
    advanced_options.add_argument('--csq_aa',
        default='Amino_acids',
        metavar='Amino_acids',
        help='Overrides the name of the CSQ field that stores amino acid change information')
    advanced_options.add_argument('--csq_pos',
        default='Protein_position',
        metavar='Protein_position',
        help='Overrides the name of the CSQ field that stores protein position information')
    advanced_options.add_argument('--csq_allele',
        default='ALLELE_NUM',
        metavar='ALLELE_NUM',
        help='Overrides the name of the CSQ field that stores allele number information')
    advanced_options.add_argument('--csq_impact',
        default='IMPACT',
        metavar='IMPACT',
        help='Overrides the name of the CSQ field that stores impact information')
    advanced_options.add_argument('--csq_class',
        default='VARIANT_CLASS',
        metavar='VARIANT_CLASS',
        help='Overrides the name of the CSQ field that stores variant class information')
    parser.add_argument('--output_basename',
        help="String to use as base name for output file")
    args = parser.parse_args()
    exit = False
    reasons = []
    if not args.protein_snvs and not args.protein_indels and not args.genomic_hotspots:
        exit = True
        reasons.append("Must provide a protein_snvs, protein_indels, or genomic_hotspots file(s)")
    if not os.path.isfile(args.vcf):
        exit = True
        reasons.append("{} is not a file.".format(args.vcf))
    if exit:
        raise FileNotFoundError("{}{}".format(chr(10),chr(10).join(reasons)))
    if not args.vcf.endswith(".vcf.gz") and not args.vcf.endswith(".vcf"):
        raise ValueError("{}{} is not a VCF or VCF.GZ".format(chr(10),args.vcf))
    else:
        return args

def csq_sanity(field,value):
    """Small sanity check for two different feature fields
    Args:
        field (str): string denoting the CSQ field
        value (str): value to compare against the field for sanity
    Return:
        str: Info about the CSQ field and hotspot value
    """
    if (field == "Gene" and not value.startswith("ENSG")) or (field == "Feature" and not value.startswith("ENST")):
        return "CSQ field: {} potentially mismatches hotspot format: {}".format(field,value)
    else:
        return "CSQ field: {}; Example hotspot: {}".format(field,value)

def process_genomic_positions(infile,hs_dict=None):
    """Function to parse a tab-delimited BED file and return a dict of its positions
    Args:
        infile (file): BED-formatted tab-delimited file containing hg38 genomic coordinates
        hs_dict (dict): dict of key: chromosome and value: list of genomic range tuples
    Return:
        dict: returns an dict list of key: chromosome and value: list of (chromStart,chromEnd) tuples
    """
    if hs_dict is None:
        hs_dict = {}
    opener = gzip.open if infile.endswith('.gz') else open
    with opener(infile,'rt') as fh:
        rd = csv.reader(fh, delimiter=chr(9))
        for row in rd:
            if not row: # skip empty rows
                continue
            if not row[0].startswith('chr'):
                continue # skip non-body rows
            chrom, chromstart, chromend = row[:3]
            if chrom not in hs_dict:
                hs_dict[chrom] = []
            if (chromstart, chromend) not in hs_dict[chrom]:
                hs_dict[chrom].append((chromstart, chromend))
    return hs_dict

def process_protein_and_position(infile,protein_colname="Hugo_Symbol",position_colname="Amino_Acid_Position",variant_class="SNV",hs_dict=None):
    """Function to parse a column-name-labeled, tab-delimited file and return protein names and positions
    Args:
        infile (file): Tab-delimited MAF file containing hotspot information
        protein_colname (string): Name of the column that contains the Hugo Symbol values
        position_colname (string): Name of the column that contains the positions of protein changes
        variant_class (string): description of the variants in the infile (SNV or INDEL)
        hs_dict (dict): dict of key: protein name and values: unique HGVSp short names
    Return:
        dict: returns key: protein name and values: unique protein positions
    """
    if hs_dict is None:
        hs_dict = {}
    opener = gzip.open if infile.endswith('.gz') else open
    with opener(infile,'rt') as fh:
        rd = csv.reader(fh, delimiter=chr(9))
        header = []
        for row in rd:
            if not row:
                continue # skip empty rows
            if row[0].startswith('#'):
                continue # skip commented rows
            if not header:
                header = row
                continue
            hugo = row[header.index(protein_colname)]
            if hugo not in hs_dict:
                hs_dict[hugo] = {"splices":[],"INDEL":[],"SNV":[]}
            pos = row[header.index(position_colname)]  # Three different types: solo positions (22), ranges (22-25), hgvsp_shorts (X128_splice)
            if 'splice' in pos: # handle splice values
                hs_dict[hugo]['splices'].append(pos)
            else:
                m = re.search(r"^(\d+|\?)-?(\d+|\?)?$", pos)
                aa_start = m[2] if m[1] == '?' else m[1]
                aa_end = m[2] if m[2] and m[2] != '?' else aa_start
                aa_range = (int(aa_start),int(aa_end)+1)
                if aa_range not in hs_dict[hugo][variant_class]:
                    hs_dict[hugo][variant_class].append(aa_range)
    return hs_dict

def create_splice_hgvsp_short(csq_record,csq_items,cons_field="Consequence",cdna_field="HGVSc"):
    """Function that will take a csq_record and the csq_items and return an HGVSp_Short value (minus the leading p.)
       if it is a splice variant and coding DNA reference sequence.
    Args:
        csq_record (str): String containing a single CSQ annotation produced by VEP
        csq_items (list): List containing the labels for the individual fields of the CSQ annotation
        cons_field (str): Name of the field in the CSQ where Consequence information can be found
        cdna_field (str): Name of the field in the CSQ where HGVSc information can be found
    Return:
        string: If the record is a coding DNA splice variant, return a modified HGVSp_Short value otherwise return an empty string
    """
    cons = csq_record.split('|')[csq_items.index(cons_field)]
    short = '' # Base value to return if we don't meet criteria
    if "splice_acceptor_variant" in cons or "splice_donor_variant" in cons: # Only perform on splice variants
        cdna_pos = csq_record.split('|')[csq_items.index(cdna_field)]
        m = re.search(r"c\.(-?\d+)",cdna_pos) # Search for coding DNA reference sequence
        if m: # If we match
            cdna_val = 1 if int(m[1]) < 0 else int(m[1]) # Rescue the negative value cDNA positions used in 5' UTRs
            pro_val = int((cdna_val + cdna_val % 3) / 3) # Recalculate protein position and drop the decimals
            short = "X{}_splice".format(pro_val)
        elif cdna_pos:
            print("Warning! Could not process HGVSc: {} of CSQ:{}{}{}".format(cdna_pos,chr(10),chr(9),csq_record))
        else:
            print("Warning! No HGVSc value in CSQ:{}{}{}".format(chr(10),chr(9),csq_record))
    return short

def process_vcf(infile,outfile,genomic_hs,hugo_hs,fid,aa='Amino_acids',pos='Protein_position',allele='ALLELE_NUM',impact='IMPACT',variant_class='VARIANT_CLASS'):
    """Function to take the input VCF and add INFO annotations for the ALT alleles that correspond to genomic or protein hotspots
    Args:
        infile (file): Path to the VCF(.GZ) file
        outfile (str): String name for the output file that will be made in the fuction
        genomic_hs (dict): Dict of hg38 genomic coordinate-based hotspots with keys: chrom and values: (chromStart, chromEnd)
        hugo_hs (dict): Dict of protein-based hotspots with keys: protein name and values: unique lists of HGVSp_short strings
        fid (str): the name of the CSQ_FIELD that is being used as the ID in the hslist; commonly Hugo Symbols: SYMBOL, Ensembl Gene IDs: Gene, or Ensembl Transcript IDs: Feature
        aa (str): the name of the CSQ_FIELD that contains the amino acid change annotation
        pos (str): the name of the CSQ_FIELD that contains the protein position annotation
        allele (str): the name of the CSQ_FIELD that contains the allele number annotation
        impact (str): the name of the CSQ_FIELD that contains the impact assessment
        variant_class (str): the name of the CSQ_FIELD that contains variant class information
    Return:
        file: A VCF.GZ file is written to Args.outfile
        int: total_hotspots_found The total number of hotspot annotations in the file
    """
    vcfin = pysam.VariantFile(infile)
    vcfin.header.info.add("HotSpotAllele",'A','Integer',"Is the corresponding ALT allele on a hotspot (0 for false; 1 for true).")
    vcfout = pysam.VariantFile(outfile,'w',header=vcfin.header)
    csq_items = vcfin.header.info["CSQ"].description.split(' Format: ')[1].split('|')
    flagged = []
    total_hotspots_found = 0
    for rec in vcfin.fetch():
        genomic_hotspot_found = False
        if rec.chrom in genomic_hs:
            for region in genomic_hs[rec.chrom]:
                region_range = range(int(region[0]),int(region[1]))
                if sufficient_overlap(range(rec.start,rec.stop),region_range,len(region_range)): # only marking exact matches
                    flagged.append("Chromosome:{} Start:{} End:{}".format(rec.chrom,region[0],region[1]))
                    genomic_hotspot_found = True
                    break
        if genomic_hotspot_found: # if it's a genomic hotspot it will be a hotspot regardless of the CSQ fields so don't check them
            total_hotspots_found += len(rec.alts)
            hotspot_array = [1] * len(rec.alts) # there should only be one alt but set for all
        else:
            hotspot_alleles = []
            for i in rec.info['CSQ']:
                field_id = i.split('|')[csq_items.index(fid)]
                imp = i.split('|')[csq_items.index(impact)]
                if imp == 'LOW': # only check if CSQ is not low impact
                    continue
                splice_short = create_splice_hgvsp_short(i,csq_items) # Get HGVSp_Short value for DNA coding splice variants
                if splice_short: # if we have a splice site we either flag it and break the loop or continue to the next CSQ without checking protein positions
                    if field_id in hugo_hs and splice_short in hugo_hs[field_id]['splices']:
                        flagged.append("Hugo_Symbol:{} HGVSp_Short:{}".format(field_id,splice_short))
                        hotspot_alleles.append(i.split('|')[csq_items.index(allele)])
                        break
                    else:
                        continue
                protein_pos = i.split('|')[csq_items.index(pos)].split('/')[0]
                if not protein_pos:
                    continue # No protein position to check. Continue to next CSQ
                elif re.search(r"^(\d+|\?)-?(\d+|\?)?$", protein_pos):
                    m = re.search(r"^(\d+|\?)-?(\d+|\?)?$", protein_pos)
                    protein_start = m[2] if m[1] == '?' else m[1]
                    protein_end = m[2] if m[2] and m[2] != '?' else protein_start
                else:
                    print("Warning! Couldn't process position: {} in CSQ:{}{}{}".format(protein_pos,chr(10),chr(9),i))
                    continue # Print anything that doesn't match our expectations for position and go to next CSQ
                protein_range = (int(protein_start),int(protein_end)+1) # correct inclusive end
                var_class = i.split('|')[csq_items.index(variant_class)]
                var_class = 'INDEL' if var_class != 'SNV' else 'SNV'
                if field_id in hugo_hs and protein_range in hugo_hs[field_id][var_class]: # Flag all exact matches
                    flagged.append("Hugo_Symbol:{} AA_Start:{} AA_End:{}".format(field_id,protein_range[0],protein_range[1]-1)) # uncorrect the inclusive addition
                    hotspot_alleles.append(i.split('|')[csq_items.index(allele)])
                    break
                elif var_class == 'INDEL' and field_id in hugo_hs: # Secondary check for indels
                    for r in hugo_hs[field_id][var_class]:
                        if protein_range[0] >= r[0] and protein_range[1] <= r[1]: # Indels just need to be contained the range
                            flagged.append("Hugo_Symbol:{} AA_Start:{} AA_End:{}".format(field_id,protein_range[0],protein_range[1]-1))
                            hotspot_alleles.append(i.split('|')[csq_items.index(allele)])
                            break
            hotspot_allele_nums = set(hotspot_alleles)
            hotspot_array = [0] * len(rec.alts)
            for i in hotspot_allele_nums:
                total_hotspots_found += 1
                hotspot_array[int(i)-1] = 1
        rec.info["HotSpotAllele"] = tuple(hotspot_array)
        vcfout.write(rec)
    if flagged:
        print("#####{}{}HotSpots Annotated:".format(chr(10),chr(10)))
        print(chr(10).join(flagged),chr(10))
    return total_hotspots_found

def build_output_name(inpath,outbase=None):
    """Builds an VCF.GZ output filename based on the input (VCF/VCF.GZ) name and given output base name
    Args:
        inpath (str): Path to the input VCF(.GZ) file
        outbase (str): Given output base name string
    Return:
        str: Filename in the format <output_basename>.<input_nameroot>.hotspots.vcf.gz
    """
    basename_split = inpath.split('/')[-1].split('.')
    inroot = basename_split[:basename_split.index('vcf')]
    if outbase:
        return outbase + '.' + '.'.join(inroot) + '.hotspots.vcf.gz'
    else:
        return '.'.join(inroot) + '.hotspots.vcf.gz'

def sufficient_overlap(x,y,overlap=2):
    """Function to determine if there is sufficient overlap between two lists
    Args:
        x (list): List containing first range
        y (list): List containing second range
        overlap (int): Minimum acceptable overlap between the two lists
    Return:
        boolean: True if there is sufficient overlap between x and y; False otherwise
    """
    return len(range(max(x[0], y[0]), min(x[-1], y[-1])+1)) >= overlap

def main():
    args = parse_command()
    genomic_hs = {}
    if args.genomic_hotspots:
        for genomic_hotspot_file in args.genomic_hotspots:
            genomic_hs = process_genomic_positions(genomic_hotspot_file,hs_dict=genomic_hs)
    hugo_hs = {}
    if args.protein_snvs:
        for hugo_hotspot_file in args.protein_snvs:
            hugo_hs = process_protein_and_position(hugo_hotspot_file,args.protein_colname,args.position_colname,variant_class="SNV",hs_dict=hugo_hs)
    if args.protein_indels:
        for hugo_hotspot_file in args.protein_indels:
            hugo_hs = process_protein_and_position(hugo_hotspot_file,args.protein_colname,args.position_colname,variant_class="INDEL",hs_dict=hugo_hs)
    output_name = build_output_name(args.vcf,args.output_basename)
    numhs = process_vcf(args.vcf,output_name,genomic_hs,hugo_hs,args.csq_field,args.csq_aa,args.csq_pos,args.csq_allele,args.csq_impact,args.csq_class)
    pysam.tabix_index(output_name, preset="vcf", force=True)
    print("Total Hotspots Annotated: {}".format(numhs))

if __name__ == '__main__':
    main()
