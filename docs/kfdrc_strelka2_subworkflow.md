# Kids First DRC Strelka2 Variant Calling Workflow
This is a subworkflow that is part of the Kids First DRC Somatic Variant Workflow that can be run as standalone.
Strelka2 is a variant caller that calls single nucleotide variants and small insertions and deletions.

![annot workflow flowchart](../docs/kfdrc_strelka2_sub_wf.cwl.png)

[This subworkflow](../sub_workflows/kfdrc_strelka2_sub_wf.cwl) does the following things as described below:

1. Run the Strelka2 variant caller tool
1. Merge the SNV and Indel results
1. Reheader merged VCF with Sample IDs provided
1. Hard filter resultant VCF on `PASS`
1. Annotate the `PASS` VCF using the [annotation sub workflow](kfdrc_annotation_subworkflow.md)
1. Rename outputs to fit a standard format

## Workflow Description and KF Recommended Inputs
This workflow runs [Strelka2 v2.9.3](https://github.com/Illumina/strelka) which calls single nucleotide variants (SNV) and insertions/deletions (INDEL).

### General Workflow Inputs
```yaml
inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  reference_dict: File
  hg38_strelka_bed: {type: File, secondaryFiles: ['.tbi']}
  manta_small_indels: {type: File?, secondaryFiles: ['.tbi']}
  use_manta_small_indels: {type: boolean?, default: false}
  input_tumor_aligned:
    type: File
    secondaryFiles: |
      ${
        var dpath = self.location.replace(self.basename, "")
        if(self.nameext == '.bam'){
          return {"location": dpath+self.nameroot+".bai", "class": "File"}
        }
        else{
          return {"location": dpath+self.basename+".crai", "class": "File"}
        }
      }
    doc: "tumor BAM or CRAM"
  input_tumor_name: string
  input_normal_aligned:
    type: File
    secondaryFiles: |
      ${
        var dpath = self.location.replace(self.basename, "")
        if(self.nameext == '.bam'){
          return {"location": dpath+self.nameroot+".bai", "class": "File"}
        }
        else{
          return {"location": dpath+self.basename+".crai", "class": "File"}
        }
      }
    doc: "normal BAM or CRAM"
  input_normal_name: string
  exome_flag: {type: ['null', string], doc: "set to 'Y' for exome mode"}
  vep_cache: {type: File, doc: "tar gzipped cache from ensembl/local converted cache"}
  vep_ref_build: {type: ['null', string], doc: "Genome ref build used, should line up with cache.", default: "GRCh38" }
  output_basename: string
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
  tool_name: {type: string?, doc: "String to describe what tool was run as part of file name", default: "strelka2_somatic"}
  # annotation vars
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots" }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots" }
  retain_info: {type: string?, doc: "csv string with INFO fields that you want to keep", default: "MQ,MQ0,QSI,HotSpotAllele"}
  retain_fmt: {type: string?, doc: "csv string with FORMAT fields that you want to keep"}
  add_common_fields: {type: boolean?, doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: true}
  bcftools_annot_columns: {type: string, doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_annot_vcf: {type: File, secondaryFiles: ['.tbi'], doc: "bgzipped annotation vcf file"}
  bcftools_public_filter: {type: string?, doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration, recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP() <= 7\"), \"AF > 0.001\"]"}
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task.", default: false }
  maf_center: {type: string?, doc: "Sequencing center of variant called", default: "."}
```
### Recommended reference inputs - all file references can be obtained [here](https://cavatica.sbgenomics.com/u/kfdrc-harmonization/kf-references/)*
Secondary files needed for each reference file will be a sub-bullet point.
* For recommendations for inputs in the `#annotation` section, see the [annotation subworkflow docs.](../sub_workflows/kfdrc_strelka2_sub_wf.cwl)
 - `indexed_reference_fasta`: `Homo_sapiens_assembly38.fasta`
   - `Homo_sapiens_assembly38.fasta.fai`
   - `Homo_sapiens_assembly38.dict`
 - `reference_dict`: `Homo_sapiens_assembly38.dict`
 - `hg38_strelka_bed`: `hg38_strelka.bed.gz`

### Source-specific inputs
 - `exome_flag`
   - `Y` if exome
   - `N` or leave blank if WGS

## Workflow outputs
```yaml
outputs:
  strelka2_prepass_vcf: {type: File, outputSource: rename_strelka_samples/reheadered_vcf}
  strelka2_protected_outputs: {type: 'File[]', outputSource: rename_protected/renamed_files}
  strelka2_public_outputs: {type: 'File[]', outputSource: rename_public/renamed_files}
```

 - `strelka2_prepass_vcf`: Combined SNV + INDEL file with renamed Sample IDs. Has all soft `FILTER` vamues generated by variatn caller
 - `strelka2_protected_outputs`: Array of files containing MAF format of PASS hits, `PASS` VCF with annotation pipeline soft `FILTER`-added values, and VCF index
 - `strelka2_public_outputs`: Same as above, except MAF and VCF have had entries with soft `FILTER` values removed

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
 | annotate_vep_annotate_vcf                                  | run step         | NA            | mkdir homo_sapiens && tar --use-compress-program="pigz -p 8" -xf /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/homo_sapiens_vep_93_GRCh38.tar.gz -C homo_sapiens && perl /ensembl-vep-release-93.7/vep --af --af_1kg --af_esp --af_gnomad --allele_number --assembly GRCh38 --biotype --buffer_size 10000 --cache --cache_version 93 --canonical --ccds --check_existing --dir_cache homo_sapiens --domains --failed 1 --fasta /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.fasta --flag_pick_allele --fork 16 --format vcf --gene_phenotype --hgvs --input_file /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_add_standard_fields/DELME_TEST.strelka2_somatic.INFO_stripped.standard.vcf.gz --no_escape --no_progress --no_stats  --numbers --offline --output_file DELME_TEST.strelka2_somatic.PASS.vep.vcf --pick_order canonical,tsl,biotype,rank,ccds,length --polyphen b --protein --pubmed  --shift_hgvs 1 --sift b --species homo_sapiens --symbol --total_length --tsl --uniprot --variant_class --vcf --xref_refseq && /ensembl-vep-release-93.7/htslib/bgzip DELME_TEST.strelka2_somatic.PASS.vep.vcf && /ensembl-vep-release-93.7/htslib/tabix DELME_TEST.strelka2_somatic.PASS.vep.vcf.gz                                                                                                                                                                                                         |
 | annotate_bcftools_gnomad_annotate                                  | run step         | NA            | bcftools annotate --annotations /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/af-only-gnomad.hg38.vcf.gz --columns INFO/AF -o DELME_TEST.strelka2_somatic.bcf_annotated.vcf.gz -O z --threads 4 /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_vep_annotate_vcf/DELME_TEST.strelka2_somatic.PASS.vep.vcf.gz && tabix DELME_TEST.strelka2_somatic.bcf_annotated.vcf.gz                                                                                                                                                                                                         |
 | annotate_gatk_add_soft_filter                                  | run step         | NA            | /gatk VariantFiltration -R /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.fasta -V /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_bcftools_gnomad_annotate/DELME_TEST.strelka2_somatic.bcf_annotated.vcf.gz -O DELME_TEST.strelka2_somatic.gatk.soft_filtered.vcf.gz --filter-name "NORM_DP_LOW" --filter-expression "vc.getGenotype('BS_30WN9M3C').getDP() <= 7" --filter-name "GNOMAD_AF_HIGH" --filter-expression "AF > 0.001"                                                                                                                                                                                                         |
 | annotate_hotspots_annotation                                  | run step         | NA            | >&2 /hotspot.py --genomic_hotspots /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/tert.bed --protein_indels /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/protein_indel_cancer_hotspots_v2.tsv --protein_snvs /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/protein_snv_cancer_hotspots_v2.tsv --output_basename DELME_TEST --vcf /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_gatk_add_soft_filter/DELME_TEST.strelka2_somatic.gatk.soft_filtered.vcf.gz                                                                                                                                                                                                         |
 | annotate_hard_filter_vcf                                  | run step         | NA            | bcftools view --include 'FILTER="PASS"|INFO/HotSpotAllele=1' /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_hotspots_annotation/DELME_TEST.DELME_TEST.strelka2_somatic.gatk.soft_filtered.hotspots.vcf.gz -O z > DELME_TEST.vcf.gz;tabix DELME_TEST.vcf.gz;                                                                                                                                                                                                         |
 | annotate_kfdrc_vcf2maf_protected                                  | run step         | NA            | gunzip -c /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_hotspots_annotation/DELME_TEST.DELME_TEST.strelka2_somatic.gatk.soft_filtered.hotspots.vcf.gz > input_file.vcf && perl /vcf2maf/vcf2maf.pl --input-vcf input_file.vcf --output-maf DELME_TEST.strelka2_somatic.vep.maf --tumor-id BS_922YMFYK --normal-id BS_30WN9M3C --ncbi-build GRCh38 --ref-fasta /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.fasta --maf-center "." --retain-info MQ,MQ0,QSI,HotSpotAllele                                                                                                                                                                                                         |
 | annotate_kfdrc_vcf2maf_public                                  | run step         | NA            | gunzip -c /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_hard_filter_vcf/DELME_TEST.vcf.gz > input_file.vcf && perl /vcf2maf/vcf2maf.pl --input-vcf input_file.vcf --output-maf DELME_TEST.strelka2_somatic.vep.maf --tumor-id BS_922YMFYK --normal-id BS_30WN9M3C --ncbi-build GRCh38 --ref-fasta /sbgenomics/Projects/0154a615-2b20-4c32-9c76-141b5eebdde1/Homo_sapiens_assembly38.fasta --maf-center "." --retain-info MQ,MQ0,QSI,HotSpotAllele                                                                                                                                                                                                         |
 | rename_protected                                  | run step         | NA            | mkdir renamed/;cp /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_hotspots_annotation/DELME_TEST.DELME_TEST.strelka2_somatic.gatk.soft_filtered.hotspots.vcf.gz renamed/DELME_TEST.strelka2_somatic.norm.annot.protected.vcf.gz;cp /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_hotspots_annotation/DELME_TEST.DELME_TEST.strelka2_somatic.gatk.soft_filtered.hotspots.vcf.gz.tbi renamed/DELME_TEST.strelka2_somatic.norm.annot.protected.vcf.gz.tbi;cp /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_kfdrc_vcf2maf_protected/DELME_TEST.strelka2_somatic.vep.maf renamed/DELME_TEST.strelka2_somatic.norm.annot.protected.maf;                                                                                                                                                                                                         |
 | rename_public                                  | run step         | NA            | mkdir renamed/;cp /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_hard_filter_vcf/DELME_TEST.vcf.gz renamed/DELME_TEST.strelka2_somatic.norm.annot.public.vcf.gz;cp /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_hard_filter_vcf/DELME_TEST.vcf.gz.tbi renamed/DELME_TEST.strelka2_somatic.norm.annot.public.vcf.gz.tbi;cp /sbgenomics/workspaces/0154a615-2b20-4c32-9c76-141b5eebdde1/tasks/731be06c-87b2-45d0-aed9-f55330240462/annotate_kfdrc_vcf2maf_public/DELME_TEST.strelka2_somatic.vep.maf renamed/DELME_TEST.strelka2_somatic.norm.annot.public.maf;                                                                                                                                                                                                         |
