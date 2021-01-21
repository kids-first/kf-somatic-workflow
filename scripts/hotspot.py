#!/usr/bin/env python

import argparse
import csv
import os
import pysam

def parse_command():
    """Function to parse the command line input

    Args:
        None

    Return:
        Namespace: returns the args as a standard Namespace object
    """
    parser = argparse.ArgumentParser(description='Flag HotSpots in VCF')
    parser.add_argument('--hotspots',
        required=True,
        help='TSV file containing the hotspots with the following format <csq-field>:<amino-acid[ref]><protein-position[start]> in the 4th column')
    parser.add_argument('--vcf',
        required=True,
        help='Input VCF or VCF.GZ file')
    parser.add_argument('--csq_field',
        default='SYMBOL',
        help='CSQ field that matches the format used in the hotspots file')
    parser.add_argument('--output_basename',
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

def process_hotspots(infile,field):
    """Function to read the hotspots file and return a list of unique hotspots

    Args:
        infile (file): file containing hotspot information
        field (str): csq field

    Return:
        list: returns a list of unique hotspot values from the file
    """
    hs_list = []
    with open(infile) as fh:
        rd = csv.reader(fh, delimiter=chr(9))
        for row in rd:
            if row[0].startswith('#'):
                continue
            hs_list.append(row[3].split(','))
    unique_hs = list(set([items for sublist in hs_list for items in sublist]))
    print(csq_sanity(field,unique_hs[0]))
    return unique_hs

def process_vcf(infile,outfile,hslist,fid,aa='Amino_acids',pos='Protein_position',allele='ALLELE_NUM'):
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
    vcfin = pysam.VariantFile(infile)
    vcfin.header.info.add("HotSpotAllele",'A','Integer',"Is the corresponding ALT allele on a hotspot (0 for false; 1 for true).")
    vcfout = pysam.VariantFile(outfile,'w',header=vcfin.header)
    csq_items = vcfin.header.info["CSQ"].description.split(' Format: ')[1].split('|')
    total_hotspots_found = 0
    for rec in vcfin.fetch():
        hotspot_allele_nums = set([i.split('|')[csq_items.index(allele)] for i in rec.info['CSQ'] if '{}:{}{}'.format(i.split('|')[csq_items.index(fid)],i.split('|')[csq_items.index(aa)].split('/')[0],i.split('|')[csq_items.index(pos)].split('/')[0]) in hslist])
### I've blown up the above line so it's easier to tell what's going on. Will keep the above for the function as I like it more
#        hotspot_alleles = []
#        for i in rec.info['CSQ']:
#            field_id = i.split('|')[csq_items.index(fid)]
#            aa_ref = i.split('|')[csq_items.index(aa)].split('/')[0]
#            pos_start = i.split('|')[csq_items.index(pos)].split('/')[0]
#            hotspot_str = '{}:{}{}'.format(field_id,aa_ref,pos_start)
#            if hotspot_str in hslist:
#                hotspot_alleles.append(i.split('|')[csq_items.index(allele))
#        hotspot_allele_nums = set(hotspot_alleles)
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

def build_output_name(inpath,outbase):
    """Builds an VCF.GZ output filename based on the input (VCF/VCF.GZ) name and given output basename

    Args:
        inpath (str): Path to the input VCF(.GZ) file
        outbase (str): Given output basename string

    Return:
        str: Filename in the format <output_basename>.<input_nameroot>.hotspots.vcf.gz
    """
    basename_split = inpath.split('/')[-1].split('.')
    inroot = basename_split[:basename_split.index('vcf')]
    return '.'.join([outbase,'.'.join(inroot) + 'hotspots.vcf.gz'])

def main():
    args = parse_command()
    uniq_hs = process_hotspots(args.hotspots,args.csq_field)
    output_name = build_output_name(args.vcf,args.output_basename)
    numhs = process_vcf(args.vcf,output_name,uniq_hs,args.csq_field)
    print("Total Hotspots Annotated: {}".format(numhs))

if __name__ == '__main__':
    main()
