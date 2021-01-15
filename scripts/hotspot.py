#!/usr/bin/env python

import argparse
import csv
import os
import pysam

"""Function to parse the command line input 

Args:
    None

Return:
    Namespace: returns the args as a standard Namespace object
"""
def parse_command():
    parser = argparse.ArgumentParser(description='Flag HotSpots in VCF')
    parser.add_argument('--hotspots',
        required=True,
        help='TSV file containing the hotspots with the following format <csq-field>:<amino-acid[ref]><protein-position[start]> in the 4th column')
    parser.add_argument('--vcf',
        required=True,
        help='Input VCF or VCF.GZ file')
    parser.add_argument('--csq_field',
        type=str,
        required=False,
        default='SYMBOL',
        help='CSQ field that matches the format used in the hotspots file')
    parser.add_argument('--output_basename',
        type=str,
        required=False,
        default="test",
        help="String to use as basename for output file")
    args = parser.parse_args()
    exit = False
    reasons = []
    if not os.path.isfile(args.hotspots):
        exit = True
        reasons.append("{} is not a file.".format(args.hotspots))
    if not os.path.isfile(args.vcf):
        exit = True
        reasons.append("{} is not a file.".format(args.vcf))
    if exit:
        raise FileNotFoundError("{}{}".format(chr(10),chr(10).join(reasons)))
    if not args.vcf.endswith(".vcf.gz") and not args.vcf.endswith(".vcf"):
        raise ValueError("{}{} is not a VCF or VCF.GZ".format(chr(10),args.vcf))
    else:
        return args

"""Small sanity check for two different feature fields

Args:
    csq (str): string denoting the CSQ field
    hotspot (str): a hotspot

Return:
    str: Info about the CSQ field and hotspot
"""
def csq_sanity(field,value):
    if (field == "Gene" and not value.startswith("ENSG")) or (field == "Feature" and not value.startswith("ENST")):
        return "CSQ field: {} potentially mismatches hotspot format: {}".format(field,value)
    else:
        return "CSQ field: {}; Example hotspot: {}".format(field,value)

"""Function to read the hotspots file and return a list of unique hotspots

Args:
    infile (file): file containing hotspot information
    csq (str): csq

Return:
    list: returns a list of unique hotspot values from the file
"""
def process_hotspots(infile,field):
    hs_list = []
    with open(infile) as fh:
        rd = csv.reader(fh, delimiter='\t')
        for row in rd:
            if row[0].startswith('#'):
                continue
            hs_list.append(row[3].split(','))
    unique_hs = list(set([items for sublist in hs_list for items in sublist]))
    print(csq_sanity(field,unique_hs[0]))
    return unique_hs

"""Function to take the input VCF and add INFO annotations for the ALT alleles that cause protein changes at hotspots

Args:
    infile (file): Path to the VCF(.GZ) file
    outfile (str): String name for the output file that will be made in the fuction
    hslist (list of strings): Unique list of hotspots derived from the provided input in the format <CSQ_FIELD>:<Amino_acids[ref]><Protein_position[start]>
    fid (str): the name of the CSQ_FIELD that is being used as the ID in the hslist; commonly Hugo Symbols: SYMBOL, Ensembl Gene IDs: Gene, or Ensembl Transcrip IDs: Feature
    aa (str): the name of the CSQ_FIELD that contains the amino acid change annotation
    pos (str): the name of the CSQ_FIELD that contains the protein position annotation
    allele (str): the name of the CSQ_FIELD that contains the allele number annotation

Return:
    file: A VCF.GZ file is written to Args.outfile
    int: total_hotspots_found The total number of hotspot annotations in the file
"""
def process_vcf(infile,outfile,hslist,fid,aa='Amino_acids',pos='Protein_position',allele='ALLELE_NUM'):
    vcfin = pysam.VariantFile(infile)
    vcfin.header.info.add("HotSpotAllele",'A','Integer',"Is the corresponding ALT allele on a hotspot (0 for false; 1 for true).")
    vcfout = pysam.VariantFile(outfile,'w',header=vcfin.header)
    csq_items = vcfin.header.info["CSQ"].description.split(' Format: ')[1].split('|')
    total_hotspots_found = 0
    for rec in vcfin.fetch():
        hotspot_allele_nums = set([i.split('|')[csq_items.index(allele)] for i in rec.info['CSQ'] if '{}:{}{}'.format(i.split('|')[csq_items.index(fid)],i.split('|')[csq_items.index(aa)].split('/')[0],i.split('|')[csq_items.index(pos)].split('/')[0]) in hslist])
        hotspot_array = [0] * len(rec.alts)
        for i in hotspot_allele_nums:
            total_hotspots_found += 1
            hotspot_array[int(i)-1] = 1
        rec.info["HotSpotAllele"] = tuple(hotspot_array)
        vcfout.write(rec)
    return total_hotspots_found

### Code block for modifying CSQ field (ultimately decided not to go this direction because CSQ is stripped immediately after in normalization)
#    for rec in vcfin.fetch():
#        new_csq = []
#        for i in rec.info['CSQ']:
#            i += "|"
#            if '{}:{}{}'.format(i.split('|')[csq_items.index(field)],i.split('|')[csq_items.index(aa)].split('/')[0],i.split('|')[csq_items.index(pos)].split('/')[0]) in hslist:
#                i += "HotSpot"
#            new_csq.append(i)
#        rec.info['CSQ'] = new_csq

"""Builds an VCF.GZ output filename based on the input (VCF/VCF.GZ) name and given output basename

Args:
    inpath (str): Path to the input VCF(.GZ) file
    outbase (srt): Given output basename string

Return:
    str: Filename in the format <output_basename>.<input_nameroot>.hotspots.vcf.gz
"""
def build_output_name(inpath,outbase):
    inroot = inpath.split('/')[-1].split('.')[:inpath.split('/')[-1].split('.').index('vcf')]
    return '.'.join([outbase,'.'.join(inroot) + 'hotspots.vcf.gz'])

def main():
    args = parse_command()
    uniq_hs = process_hotspots(args.hotspots,args.csq_field)
    output_name = build_output_name(args.vcf,args.output_basename)
    numhs = process_vcf(args.vcf,output_name,uniq_hs,args.csq_field)
    print("Total Hotspots Annotated: {}".format(numhs))

if __name__ == '__main__':
    main()