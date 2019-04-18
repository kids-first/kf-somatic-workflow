# KFDRC Somatic Caller Workflow
Work in progress somatic calling workflow

## Introduction
Currently produces somatic snv and indel calls from strelka2  and mutect2, structural variation calls from manta, and copy number variation from ControlFreec.  GATK4 SelectVariants used to filter vcfs on `PASS`. Variant effect predictor and MSKCC's vcf2maf used to annotate output vcfs and generate maf files.

### Strelka2
v2.9.3, somatic snv and indel calls, https://github.com/Illumina/strelka
### Manta
v1.4.0, structural varation calls, https://github.com/Illumina/manta
### Control-FREEC
v11.5, copynumber variation https://github.com/BoevaLab/FREEC
### GATK
v4.1.1, used to run Mutect2 and perform common tasks, like merge vcfs from callers.  Mutect2 running and filtering follows strictly Broad workflow form https://github.com/broadinstitute/gatk/blob/4.1.1.0/scripts/mutect2_wdl/mutect2.wdl
### Variant effect predictor
v93, used to annotate vcfs and generate mafs per MSKCC specs https://github.com/Ensembl/ensembl-vep
### vcf2maf
v1.6.16 https://github.com/mskcc/vcf2maf/releases


## Usage

### Inputs:
```yaml
inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  reference_dict: File
  wgs_calling_interval_list: File
  af_only_gnomad_vcf: {type: File, secondaryFiles: ['.tbi']}
  exac_common_vcf: {type: File, secondaryFiles: ['.tbi']}
  hg38_strelka_bed: File
  input_tumor_aligned: File
  input_tumor_name: string
  input_normal_aligned: File
  input_normal_name: string
  threads: {type: int, doc: "For ControlFreeC.  Recommend 16 max, as I/O gets saturated after that losing any advantage."}
  exome_flag: ['null', string]
  select_vars_mode: {type: string, doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression"}
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  chr_len: File
  ref_chrs: File
  output_basename: string
```

### Suggested inputs:
```text
  indexed_reference_fasta: Homo_sapiens_assembly38.fasta
  reference_dict: Homo_sapiens_assembly38.dict
  wgs_calling_interval_list: wgs_canonical_calling_regions.hg38.interval_list
  af_only_gnomad_vcf: af-only-gnomad.hg38.vcf.gz
  exac_common_vcf: small_exac_common_3.hg38.vcf.gz
  hg38_strelka_bed: hg38_strelka.bed.gz
  threads: 16
  select_vars_mode: gatk
  vep_cache: homo_sapiens_vep_93_GRCh38_convert_cache.tar.gz
  chr_len: hs38_chr.len
  ref_chrs: GRCh38_everyChrs.tar.gz
  ```

  ### Outputs:
  ```yaml
  outputs:
    strelka2_vep_vcf: {type: File, outputSource: vep_annot_strelka2/output_vcf}
    strelka2_vep_tbi: {type: File, outputSource: vep_annot_strelka2/output_tbi}
    strelka2_prepass_vcf: {type: File, outputSource: rename_strelka_samples/reheadered_vcf}
    strelka2_vep_maf: {type: File, outputSource: vep_annot_strelka2/output_maf}
    mutect2_vep_vcf: {type: File, outputSource: vep_annot_mutect2/output_vcf}
    mutect2_vep_tbi: {type: File, outputSource: vep_annot_mutect2/output_tbi}
    mutect2_prepass_vcf: {type: File, outputSource: filter_mutect2_vcf/filtered_vcf}
    mutect2_vep_maf: {type: File, outputSource: vep_annot_mutect2/output_maf}
    manta_vep_vcf: {type: File, outputSource: vep_annot_manta/output_vcf}
    manta_vep_tbi: {type: File, outputSource: vep_annot_manta/output_tbi}
    manta_prepass_vcf: {type: File, outputSource: rename_manta_samples/reheadered_vcf}
    manta_vep_maf: {type: File, outputSource: vep_annot_manta/output_maf}
    cnv_bam_ratio: { type: File, outputSource: control_free_c/output_txt }
    cnv_pval: { type: File, outputSource: control_free_c_r/output_pval }
    cnv_png: { type: File, outputSource: control_free_c_viz/output_png }
```
