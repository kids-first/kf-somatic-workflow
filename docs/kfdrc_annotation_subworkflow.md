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

## Workflow outputs
```yaml
outputs:
  annotated_protected_vcf: {type: File, outputSource: hotspots_annotation/hotspots_vcf}
  annotated_protected_maf: {type: File, outputSource: kfdrc_vcf2maf_protected/output_maf}
  annotated_public_vcf: {type: File, outputSource: hard_filter_vcf/filtered_vcf}
  annotated_public_maf: {type: File, outputSource: kfdrc_vcf2maf_public/output_maf}
```
 - `annotated_protected_vcf`: `PASS` VCF with annotation pipeline soft `FILTER`-added values, and VCF index
 - `annotated_protected_maf`: MAF format of `annotated_protected_vcf`
 - `annotated_public_vcf`: Same as `annotated_protected_vcf`, hard-filtered to include `PASS` only
 - `annotated_public_maf`: MAF format of `annotated_public_vcf`

# Simulated bash calls
For your convenience, we have exported the calls from an example workflow, as they were run in each docker image for each step:


 | Step                         | Type         | Num scatter            | Command                                                                                                                                                                                                         |
 | ---------------------------- | ------------ | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
 | normalize_vcf                         | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | normalize_vcf                         | run step         | NA            | VCF=/sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/DELME_ANNOT_TEST/e8ae0f0f-20dd-447c-b04c-63b53d867d7d.mutect2_somatic.PASS.vep.vcf.gz                                                                                                                                                                                                         |
 | normalize_vcf                         | run step         | NA            |  >&2 echo checking if strip flag given; >&2 echo no strip flag given && bcftools norm -m '-any' $VCF > DELME_ANNOT_TEST.mutect2_somatic.bcf_norm.vcf && /vt/vt normalize -n -r /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.fasta DELME_ANNOT_TEST.mutect2_somatic.bcf_norm.vcf > DELME_ANNOT_TEST.mutect2_somatic.bcf_vt_norm.vcf && bgzip DELME_ANNOT_TEST.mutect2_somatic.bcf_vt_norm.vcf && tabix DELME_ANNOT_TEST.mutect2_somatic.bcf_vt_norm.vcf.gz                                                                                                                                                                                                         |
 | bcftools_strip_info                         | run step         | NA            | /bin/bash -c set -eo pipefail                                                                                                                                                                                                         |
 | bcftools_strip_info                         | run step         | NA            | (bcftools annotate -x INFO/AF /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/584d7621-658c-46f8-b0cd-7ac49e788549/normalize_vcf/DELME_ANNOT_TEST.mutect2_somatic.bcf_vt_norm.vcf.gz -O z  -o DELME_ANNOT_TEST.mutect2_somatic.INFO_stripped.vcf.gz && tabix DELME_ANNOT_TEST.mutect2_somatic.INFO_stripped.vcf.gz) || (echo "Check errors, likely does not have INFO, trying to pass input instead" >&2; cp /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/584d7621-658c-46f8-b0cd-7ac49e788549/normalize_vcf/DELME_ANNOT_TEST.mutect2_somatic.bcf_vt_norm.vcf.gz .; cp /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/584d7621-658c-46f8-b0cd-7ac49e788549/normalize_vcf/DELME_ANNOT_TEST.mutect2_somatic.bcf_vt_norm.vcf.gz.tbi  .;)                                                                                                                                                                                                         |
 | add_standard_fields                         | run step         | NA            | >&2 echo 'User opted to skip adding fields to vcf' && exit 0; --strelka2_vcf /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/584d7621-658c-46f8-b0cd-7ac49e788549/bcftools_strip_info/DELME_ANNOT_TEST.mutect2_somatic.INFO_stripped.vcf.gz --tumor_name BS_79SYEHY3 --normal_name BS_0KKH9VKP --output_basename DELME_ANNOT_TEST                                                                                                                                                                                                         |
 | vep_annotate_vcf                         | run step         | NA            | mkdir homo_sapiens && tar --use-compress-program="pigz -p 8" -xf /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/homo_sapiens_vep_93_GRCh38.tar.gz -C homo_sapiens && perl /ensembl-vep-release-93.7/vep --af --af_1kg --af_esp --af_gnomad --allele_number --assembly GRCh38 --biotype --buffer_size 10000 --cache --cache_version 93 --canonical --ccds --check_existing --dir_cache homo_sapiens --domains --failed 1 --fasta /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.fasta --flag_pick_allele --fork 16 --format vcf --gene_phenotype --hgvs --input_file /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/584d7621-658c-46f8-b0cd-7ac49e788549/add_standard_fields/DELME_ANNOT_TEST.mutect2_somatic.INFO_stripped.vcf.gz --no_escape --no_progress --no_stats  --numbers --offline --output_file DELME_ANNOT_TEST.mutect2_somatic.PASS.vep.vcf --pick_order canonical,tsl,biotype,rank,ccds,length --polyphen b --protein --pubmed  --shift_hgvs 1 --sift b --species homo_sapiens --symbol --total_length --tsl --uniprot --variant_class --vcf --xref_refseq && /ensembl-vep-release-93.7/htslib/bgzip DELME_ANNOT_TEST.mutect2_somatic.PASS.vep.vcf && /ensembl-vep-release-93.7/htslib/tabix DELME_ANNOT_TEST.mutect2_somatic.PASS.vep.vcf.gz                                                                                                                                                                                                         |
 | bcftools_gnomad_annotate                         | run step         | NA            | bcftools annotate --annotations /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/af-only-gnomad.hg38.vcf.gz --columns INFO/AF -o DELME_ANNOT_TEST.mutect2_somatic.bcf_annotated.vcf.gz -O z --threads 4 /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/584d7621-658c-46f8-b0cd-7ac49e788549/vep_annotate_vcf/DELME_ANNOT_TEST.mutect2_somatic.PASS.vep.vcf.gz && tabix DELME_ANNOT_TEST.mutect2_somatic.bcf_annotated.vcf.gz                                                                                                                                                                                                         |
 | gatk_add_soft_filter                         | run step         | NA            | /gatk VariantFiltration -R /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.fasta -V /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/584d7621-658c-46f8-b0cd-7ac49e788549/bcftools_gnomad_annotate/DELME_ANNOT_TEST.mutect2_somatic.bcf_annotated.vcf.gz -O DELME_ANNOT_TEST.mutect2_somatic.gatk.soft_filtered.vcf.gz --filter-name "NORM_DP_LOW" --filter-expression "vc.getGenotype('BS_0KKH9VKP').getDP() <= 7" --filter-name "GNOMAD_AF_HIGH" --filter-expression "AF > 0.001"                                                                                                                                                                                                         |
 | hotspots_annotation                         | run step         | NA            | >&2 /hotspot.py --genomic_hotspots /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/tert.bed --protein_indels /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/protein_indel_cancer_hotspots_v2.tsv --protein_snvs /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/protein_snv_cancer_hotspots_v2.tsv --output_basename DELME_ANNOT_TEST --vcf /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/584d7621-658c-46f8-b0cd-7ac49e788549/gatk_add_soft_filter/DELME_ANNOT_TEST.mutect2_somatic.gatk.soft_filtered.vcf.gz                                                                                                                                                                                                         |
 | hard_filter_vcf                         | run step         | NA            | bcftools view --include 'FILTER="PASS"|INFO/HotSpotAllele=1' /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/584d7621-658c-46f8-b0cd-7ac49e788549/hotspots_annotation/DELME_ANNOT_TEST.DELME_ANNOT_TEST.mutect2_somatic.gatk.soft_filtered.hotspots.vcf.gz -O z > DELME_ANNOT_TEST.vcf.gz;tabix DELME_ANNOT_TEST.vcf.gz;                                                                                                                                                                                                         |
 | kfdrc_vcf2maf_protected                         | run step         | NA            | gunzip -c /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/584d7621-658c-46f8-b0cd-7ac49e788549/hotspots_annotation/DELME_ANNOT_TEST.DELME_ANNOT_TEST.mutect2_somatic.gatk.soft_filtered.hotspots.vcf.gz > input_file.vcf && perl /vcf2maf/vcf2maf.pl --input-vcf input_file.vcf --output-maf DELME_ANNOT_TEST.mutect2_somatic.vep.maf --tumor-id BS_79SYEHY3 --normal-id BS_0KKH9VKP --ncbi-build GRCh38 --ref-fasta /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.fasta --maf-center "." --retain-info MBQ,TLOD,HotSpotAllele                                                                                                                                                                                                         |
 | kfdrc_vcf2maf_public                         | run step         | NA            | gunzip -c /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/584d7621-658c-46f8-b0cd-7ac49e788549/hard_filter_vcf/DELME_ANNOT_TEST.vcf.gz > input_file.vcf && perl /vcf2maf/vcf2maf.pl --input-vcf input_file.vcf --output-maf DELME_ANNOT_TEST.mutect2_somatic.vep.maf --tumor-id BS_79SYEHY3 --normal-id BS_0KKH9VKP --ncbi-build GRCh38 --ref-fasta /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.fasta --maf-center "." --retain-info MBQ,TLOD,HotSpotAllele                                                                                                                                                                                                         |
