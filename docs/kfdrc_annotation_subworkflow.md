# Kids First DRC Somatic Variant Annotation Workflow
This is a subworkflow that is part of the Kids First DRC Somatic Variant Workflow that can be run as standalone.
Annotation of variant calls helps give context to the possible biological consequences of each variant.

![annot workflow flowchart](../docs/somatic_annotation_wf.png)

It does the following things as described below:

1. Normalize VCF
1. Strip specified `INFO` and `FORMAT` fields (Only if adding a new annotation that clashes with existing)
1. Annotate with VEP
1. Annotated with an additional vcf - optional, recommend using AF-only WGS gnomad VCF)
1. Soft filter on remarkable variant characteristics
   - KF recommends normal read depth <= 7 and gnomAD AF > 0.001
   - This output will be considered `protected`
1. Annotate with hotspots - KF recommends cancer genome hotspots v2, formatting required and explained below
1. Create MAF output using a modified version of MSKCC's vcf2maf
1. Hard filter on vcf based on user-specified criteria - this output would be considered `public`
1. Create MAF output based on `public` vcf

## Workflow Description and KF Recommended Inputs
The additional gnomAD annotation, hotspot annotation, and soft + hard filtering are part of process called "Germline Masking."
The purpose of this is to create outputs that are safe for public consumption by marking and creating a version of outputs deemed a "germline risk" based on specified criteria.
For KF, based on guidance from the Genomic Data Commons (GDC), this means filtering variants with a normal read depth of <= 7 reads and a gnomAD AF > 0.001.
The gnomAD AF filter is pretty intuitive - gnomAD is a database resource of variants and their estimated prevalence in the human population.
Therefore, a variant that is higher than the recommended threshold can be seen as a higher risk of being a common and identifiable variant, and a lower risk for being disease-causing.
The normal depth argument may be less intuitive to some, and an example might help explain its importance:
You've got a somatic variant candidate, and the normal-sample reads all support the reference allele.
However,  if there are only ~5 of those reads, then there's a decent chance that it's a germline het variant and you just didn't get lucky enough to see any alt reads in the normal sample.
A prior understanding that heterozygous germline variant are much more common than somatic variants informs this.

### General workflow inputs:
```yaml
inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  input_vcf: {type: File, secondaryFiles: ['.tbi'], doc: "Input vcf to annotate and soft filter"}
  input_tumor_name: string
  input_normal_name: string
  add_common_fields: {type: boolean, doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: false}
  bcftools_annot_columns: {type: string, doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF"}
  bcftools_annot_vcf: {type: File, secondaryFiles: ['.tbi'], doc: "bgzipped annotation vcf file"}
  bcftools_public_filter: {type: string?, doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration for clues"}
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  vep_ref_build: {type: ['null', string], doc: "Genome ref build used, should line up with cache.", default: "GRCh38" }
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task." }
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots" }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots" }
  output_basename: string
  tool_name: string
  retain_info: {type: string?, doc: "csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,Hotspot`"}
  retain_fmt: {type: string?, doc: "csv string with FORMAT fields that you want to keep"}
  maf_center: {type: string?, doc: "Sequencing center of variant called", default: "."}
```

### Recommended reference inputs - all file references can be obtained [here](https://cavatica.sbgenomics.com/u/kfdrc-harmonization/kf-references/)
Secondary files needed for each reference file will be a sub-bullet point
 - `indexed_reference_fasta`: `Homo_sapiens_assembly38.fasta`
   - `Homo_sapiens_assembly38.fasta.fai`
   - `Homo_sapiens_assembly38.dict`
 - `bcftools_annot_columns`: "INFO/AF"
 - `bcftools_annot_vcf`: `af-only-gnomad.hg38.vcf.gz`
   - `af-only-gnomad.hg38.vcf.gz.tbi`
 - `bcftools_public_filter`: 'FILTER="PASS"|INFO/HotSpotAllele=1'
 - `gatk_filter_name`: ["NORM_DP_LOW", "GNOMAD_AF_HIGH"]
 - `gatk_filter_expression`: ["vc.getGenotype('insert_normal_sample_name').getDP() <= 7", "AF > 0.001"] # NOTE!! Replace `insert_normal_sample_name` with the value you'd use for `input_normal_name`!
 - `vep_cache`: `homo_sapiens_vep_93_GRCh38.tar.gz`
 - `genomic_hotspots`: `tert.bed` # This file has two common TERT promoter gene hot spots
 - `protein_snv_hotspots`: `protein_snv_cancer_hotspots_v2.tsv` # A tsv formatted SNV + MNV subset of https://www.cancerhotspots.org/files/hotspots_v2.xls
 - `protein_indel_hotspots`: `protein_indel_cancer_hotspots_v2.tsv` # A tsv formatted INDEL subset of https://www.cancerhotspots.org/files/hotspots_v2.xls

### Source-specific inputs
For each input, the sub-bullet refers to when to use the suggested input
 - `add_common_fields`
   - Strelka2 calls: `true`
   - All others: `false`
 - `retain_info` # This is fairly subjective, some useful columns unique from each caller to carry over from VCF to MAF
   - Strelka2: "MQ,MQ0,QSI,HotSpotAllele"
   - Mutect2: "MBQ,TLOD,HotSpotAllele"
   - Lancet: "MS,FETS,HotSpotAllele"
   - Vardict: "MSI,MSILEN,SOR,SSF,HotSpotAllele"
 - `tool_name`:
   - `Strelka2`: `strelka2_somatic`
   - `Mutect2`: `mutect2_somatic`
   - `Lancet`: `lancet_somatic`
   - `VarDict Java`: `vardict_somatic`

# Simulated Bash Script
As a convenience, we have extracted the commands run within each docker image during the workflow

 | Step                                  | Type         | Num scatter            | Command                                                                                                                                                                                                         |
 | ------------------------------------- | ------------ | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
 | strelka2                                  | run step         | NA            | /strelka-2.9.3.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/NORMAL_SUBSET_TEST.cram --tumorBam /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/TUMOR_SUBSET_TEST.cram --ref /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.fasta --callRegions /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/hg38_strelka.bed.gz  --runDir=./ && ./runWorkflow.py -m local -j 18                                                                                                                                                                                                         |
 | merge_strelka2_vcf                                  | run step         | NA            | /gatk MergeVcfs --java-options "-Xmx2000m" --TMP_DIR=./TMP --CREATE_INDEX=true --SEQUENCE_DICTIONARY=/sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.dict --OUTPUT=DELME_TEST.strelka2_somatic.merged.vcf.gz  -I /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/strelka2/results/variants/somatic.snvs.vcf.gz -I /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/strelka2/results/variants/somatic.indels.vcf.gz                                                                                                                                                                                                         |
 | rename_strelka_samples                                  | run step         | NA            | echo BS_30WN9M3C > sample_list.txt && echo BS_922YMFYK >> sample_list.txt && bcftools reheader -s sample_list.txt /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/merge_strelka2_vcf/DELME_TEST.strelka2_somatic.merged.vcf.gz > DELME_TEST.strelka2_somatic.merged.reheadered.vcf.gz && tabix DELME_TEST.strelka2_somatic.merged.reheadered.vcf.gz                                                                                                                                                                                                         |
 | gatk_selectvariants_strelka2                                  | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | gatk_selectvariants_strelka2                                  | run step         | NA            | /gatk SelectVariants --java-options "-Xmx8000m" -V /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/rename_strelka_samples/DELME_TEST.strelka2_somatic.merged.reheadered.vcf.gz -O DELME_TEST.strelka2_somatic.PASS.vcf.gz --exclude-filtered TRUE                                                                                                                                                                                                         |
 | annotate_normalize_vcf                                  | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | annotate_normalize_vcf                                  | run step         | NA            | VCF=/sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/gatk_selectvariants_strelka2/DELME_TEST.strelka2_somatic.PASS.vcf.gz                                                                                                                                                                                                         |
 | annotate_normalize_vcf                                  | run step         | NA            |  >&2 echo checking if strip flag given; >&2 echo no strip flag given && bcftools norm -m '-any' $VCF > DELME_TEST.strelka2_somatic.bcf_norm.vcf && /vt/vt normalize -n -r /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.fasta DELME_TEST.strelka2_somatic.bcf_norm.vcf > DELME_TEST.strelka2_somatic.bcf_vt_norm.vcf && bgzip DELME_TEST.strelka2_somatic.bcf_vt_norm.vcf && tabix DELME_TEST.strelka2_somatic.bcf_vt_norm.vcf.gz                                                                                                                                                                                                         |
 | annotate_bcftools_strip_info                                  | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | annotate_bcftools_strip_info                                  | run step         | NA            | (bcftools annotate -x INFO/AF /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_normalize_vcf/DELME_TEST.strelka2_somatic.bcf_vt_norm.vcf.gz -O z  -o DELME_TEST.strelka2_somatic.INFO_stripped.vcf.gz && tabix DELME_TEST.strelka2_somatic.INFO_stripped.vcf.gz) || (echo "Check errors, likely does not have INFO, trying to pass input instead" >&2; cp /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_normalize_vcf/DELME_TEST.strelka2_somatic.bcf_vt_norm.vcf.gz .; cp /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_normalize_vcf/DELME_TEST.strelka2_somatic.bcf_vt_norm.vcf.gz.tbi  .;)                                                                                                                                                                                                         |
 | annotate_add_standard_fields                                  | run step         | NA            | >&2 /usr/bin/add_strelka2_fields.py --strelka2_vcf /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_bcftools_strip_info/DELME_TEST.strelka2_somatic.INFO_stripped.vcf.gz --tumor_name BS_922YMFYK --normal_name BS_30WN9M3C --output_basename DELME_TEST                                                                                                                                                                                                         |
 | annotate_vep_annotate_vcf                                  | run step         | NA            | mkdir homo_sapiens && tar --use-compress-program="pigz -p 8" -xf /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/homo_sapiens_vep_93_GRCh38.tar.gz -C homo_sapiens && perl /ensembl-vep-release-93.7/vep --af --af_1kg --af_esp --af_gnomad --allele_number --assembly GRCh38 --biotype --buffer_size 10000 --cache --cache_version 93 --canonical --ccds --check_existing --dir_cache homo_sapiens --domains --failed 1 --fasta /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.fasta --flag_pick_allele --fork 16 --format vcf --gene_phenotype --hgvs --input_file /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_add_standard_fields/DELME_TEST.strelka2_somatic.INFO_stripped.standard.vcf.gz --no_escape --no_progress --no_stats  --numbers --offline --output_file DELME_TEST.strelka2_somatic.PASS.vep.vcf --pick_order canonical,tsl,biotype,rank,ccds,length --polyphen b --protein --pubmed  --shift_hgvs 1 --sift b --species homo_sapiens --symbol --total_length --tsl --uniprot --variant_class --vcf --xref_refseq && /ensembl-vep-release-93.7/htslib/bgzip DELME_TEST.strelka2_somatic.PASS.vep.vcf && /ensembl-vep-release-93.7/htslib/tabix DELME_TEST.strelka2_somatic.PASS.vep.vcf.gz                                                                                                                                                                                                        