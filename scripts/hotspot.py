#!/usr/bin/env python3

import argparse
import csv
import os
import pysam
import gzip

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
        help="Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots")
    hotspot_options.add_argument('--protein_hotspots',
        nargs='*',
        help="Column-name-containing, tab-delimited file(s) containing protein names and HGVSp short values corresponding to hotspots")
    required_options.add_argument('--vcf',
        required=True,
        help='Input VCF or VCF.GZ file')
    advanced_options.add_argument('--protein_colname',
        default='Hugo_Symbol',
        metavar='Hugo_Symbol',
        help='Overrides the column name in the protein_hotspots file(s) that contains the protein name information')
    advanced_options.add_argument('--hgvsp_colname',
        default='HGVSp_Short',
        metavar='HGVSp_Short',
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
    parser.add_argument('--output_basename',
        help="String to use as basename for output file")
    args = parser.parse_args()
    exit = False
    reasons = []
    if not args.protein_hotspots and not args.genomic_hotspots:
        exit = True
        reasons.append("Must provide either protein_hotspots or genomic_hotspots file(s)")
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

def process_protein_and_hgvsp(infile,protein_colname="Hugo_Symbol",hgvsp_colname="HGVSp_Short",hs_dict=None):
    """Function to parse a column-name-labeled, tab-delimited file and return protein name and HGVSp
    Args:
        infile (file): Tab-delimited MAF file containing hotspot information
        protein_colname (string): Name of the column that contains the Hugo Symbol values
        hgvsp_colname (string): Name of the column that contains the HGVSp Short values
        hs_dict (dict): dict of key: protein name and values: unique HGVSp short names
    Return:
        dict: returns key: protein name and values: unique HGVSp short names
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
            hgvsp = row[header.index(hgvsp_colname)]  #('A2ML1', 'p.V598=')
            if hugo not in hs_dict:
                hs_dict[hugo] = []
            if hgvsp not in hs_dict[hugo]:
                hs_dict[hugo].append(hgvsp)
    return hs_dict

def process_vcf(infile,outfile,genomic_hs,hugo_hs,fid,aa='Amino_acids',pos='Protein_position',allele='ALLELE_NUM',impact='IMPACT'):
    """Function to take the input VCF and add INFO annotations for the ALT alleles that correspond to genomic or protein hotspots
    Args:
        infile (file): Path to the VCF(.GZ) file
        outfile (str): String name for the output file that will be made in the fuction
        genomic_hs (dict): Dict of hg38 genomic coordinate-based hotspots with keys: chrom and values: (chromStart, chromEnd)
        hugo_hs (dict): Dict of protein-based hotspots with keys: protein name and values: unique lists of HGVSp_short strings
        fid (str): the name of the CSQ_FIELD that is being used as the ID in the hslist; commonly Hugo Symbols: SYMBOL, Ensembl Gene IDs: Gene, or Ensembl Transcrip IDs: Feature
        aa (str): the name of the CSQ_FIELD that contains the amino acid change annotation
        pos (str): the name of the CSQ_FIELD that contains the protein position annotation
        allele (str): the name of the CSQ_FIELD that contains the allele number annotation
        impact (str): the name of the CSQ_FIELD that contains the impact assessment
    Return:
        file: A VCF.GZ file is written to Args.outfile
        int: total_hotspots_found The total number of hotspot annotations in the file
    """
    vcfin = pysam.VariantFile(infile)
    vcfin.header.info.add("HotSpotAllele",'A','Integer',"Is the corresponding ALT allele on a hotspot (0 for false; 1 for true).")
    vcfout = pysam.VariantFile(outfile,'w',header=vcfin.header)
    csq_items = vcfin.header.info["CSQ"].description.split(' Format: ')[1].split('|')
    total_hotspots_found = 0
    for rec in vcfin.fetch():
        genomic_hotspot_found = False
        if rec.chrom in genomic_hs:
            for region in genomic_hs[rec.chrom]:
                if rec.start == int(region[0]) and rec.stop == int(region[1]): # we're only marking it if the rec start and stop match exactly to the genomic regions
                    genomic_hotspot_found = True
                    break
        if genomic_hotspot_found: # if it's a genomic hotspot it will be a hotspot regardless of the CSQ fields so don't check them
            total_hotspots_found += len(rec.alts)
            hotspot_array = [1] * len(rec.alts) # there should only be one alt but set for all
        else:
            hotspot_alleles = []
            for i in rec.info['CSQ']:
                field_id = i.split('|')[csq_items.index(fid)]
                aa_ref = i.split('|')[csq_items.index(aa)].split('/')[0]
                pos_start = i.split('|')[csq_items.index(pos)].split('/')[0]
                imp = i.split('|')[csq_items.index(impact)]
                if imp is 'LOW':
                    continue # If the impact is low we won't flag it so skip the HUGO checks
                aa_pos = aa_ref + pos_start
                if field_id in hugo_hs:
                    for hgvsp in hugo_hs[field_id]:
                        if aa_pos in hgvsp: # check the aa_pos partial string against all the hgvsps string for that protein
                            hotspot_alleles.append(i.split('|')[csq_items.index(allele)])
                            break # that CSQ has a hotspot so you can stop checking hugos and move on to the next vcf rec
            hotspot_allele_nums = set(hotspot_alleles)
            hotspot_array = [0] * len(rec.alts)
            for i in hotspot_allele_nums:
                total_hotspots_found += 1
                hotspot_array[int(i)-1] = 1
        rec.info["HotSpotAllele"] = tuple(hotspot_array)
        vcfout.write(rec)
    return total_hotspots_found

def build_output_name(inpath,outbase=None):
    """Builds an VCF.GZ output filename based on the input (VCF/VCF.GZ) name and given output basename
    Args:
        inpath (str): Path to the input VCF(.GZ) file
        outbase (str): Given output basename string
    Return:
        str: Filename in the format <output_basename>.<input_nameroot>.hotspots.vcf.gz
    """
    basename_split = inpath.split('/')[-1].split('.')
    inroot = basename_split[:basename_split.index('vcf')]
    if outbase:
        return outbase + '.' + '.'.join(inroot) + '.hotspots.vcf.gz'
    else:
        return '.'.join(inroot) + '.hotspots.vcf.gz'

def main():
    args = parse_command()
    genomic_hs = {}
    if args.genomic_hotspots:
        for genomic_hotspot_file in args.genomic_hotspots:
            genomic_hs = process_genomic_positions(genomic_hotspot_file,hs_dict=genomic_hs)
    hugo_hs = {}
    if args.protein_hotspots:
        for hugo_hotspot_file in args.protein_hotspots:
            hugo_hs = process_protein_and_hgvsp(hugo_hotspot_file,args.protein_colname,args.hgvsp_colname,hs_dict=hugo_hs)
    output_name = build_output_name(args.vcf,args.output_basename)
    numhs = process_vcf(args.vcf,output_name,genomic_hs,hugo_hs,args.csq_field,args.csq_aa,args.csq_pos,args.csq_allele,args.csq_impact)
    pysam.tabix_index(output_name, preset="vcf", force=True) 
    print("Total Hotspots Annotated: {}".format(numhs))

if __name__ == '__main__':
    main()
