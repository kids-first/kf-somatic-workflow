cwlVersion: v1.0
class: Workflow
id: kfdrc_production_strelka2_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  # Required
  reference_fasta: { type: File }
  reference_fai: { type: 'File?' }
  reference_dict: { type: 'File?' }
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
  manta_small_indels: { type: 'File?', doc: "Small indels output from Manta workflow", secondaryFiles: [.tbi] }
  vep_cache: { type: File, doc: "tar gzipped cache from ensembl/local converted cache" }
  hg38_strelka_bed: { type: File, doc: "Bgzipped interval bed file. Recommned padding 100bp for WXS; Recommend canonical chromosomes for WGS" }
  hg38_strelka_tbi: { type: 'File?', doc: "Tabix index for hg38_strelka_bed" }
  output_basename: { type: string, doc: "String value to use as basename for outputs" }
  wgs_or_wxs: { type: { type: enum, name: sex, symbols: ["WGS", "WXS"] }, doc: "Select if this run is WGS or WXS" }

  # Optional with One Default
  select_vars_mode: { type: ['null', { type: enum, name: select_vars_mode, symbols: ["gatk", "grep"] }], default: "gatk", doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression" }
  vep_ref_build: { type: 'string?', default: "GRCh38", doc: "Genome ref build used, should line up with cache" }

  # Optional with Multiple Defaults (handled in choose_defaults)
  exome_flag: { type: 'string?', doc: "Whether to run in exome mode for callers. Y for WXS, N for WGS" }

outputs:
  strelka2_vep_vcf: { type: File, outputSource: run_strelka2/strelka2_vep_vcf }
  strelka2_vep_tbi: { type: File, outputSource: run_strelka2/strelka2_vep_tbi }
  strelka2_prepass_vcf: { type: File, outputSource: run_strelka2/strelka2_prepass_vcf }
  strelka2_vep_maf: { type: File, outputSource: run_strelka2/strelka2_vep_maf }

steps:
  choose_defaults:
    run: ../tools/mode_defaults.cwl
    in:
      input_mode: wgs_or_wxs
      exome_flag: exome_flag
    out: [out_exome_flag,out_cnvkit_wgs_mode,out_i_flag,out_lancet_padding,out_lancet_window,out_vardict_padding]

  prepare_reference:
    run: ../sub_workflows/prepare_reference.cwl
    in:
      input_fasta: reference_fasta
      input_fai: reference_fai
      input_dict: reference_dict
    out: [indexed_fasta,reference_dict]

  index_strelka_bed:
    run: ../tools/tabix_index.cwl
    in:
      input_file: hg38_strelka_bed
      input_index: hg38_strelka_tbi
    out: [output]

  run_strelka2:
    run: ../sub_workflows/kfdrc_strelka2_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      reference_dict: prepare_reference/reference_dict
      hg38_strelka_bed: index_strelka_bed/output
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      manta_small_indels: manta_small_indels
      exome_flag: exome_flag
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      output_basename: output_basename
      select_vars_mode: select_vars_mode
    out:
      [strelka2_vep_vcf, strelka2_vep_tbi, strelka2_prepass_vcf, strelka2_vep_maf]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 6
  - class: 'sbg:AWSInstanceType'
    value: c5.9xlarge
