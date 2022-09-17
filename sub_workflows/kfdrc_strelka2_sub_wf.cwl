cwlVersion: v1.0
class: Workflow
id: kfdrc_strelka2_sub_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  indexed_reference_fasta: {type: 'File', secondaryFiles: [.fai, ^.dict]}
  reference_dict: File
  hg38_strelka_bed: {type: 'File', secondaryFiles: ['.tbi']}
  manta_small_indels: {type: 'File?', secondaryFiles: ['.tbi']}
  use_manta_small_indels: {type: 'boolean?', default: false}
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
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
  tool_name: {type: 'string?', doc: "String to describe what tool was run as part of file name", default: "strelka2_somatic"}
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache"}
  vep_ram: {type: 'int?', doc: "In GB, may need to increase this value depending on the size/complexity of input"}
  vep_cores: {type: 'int?', doc: "Number of cores to use. May need to increase for really large inputs"}
  vep_buffer_size: {type: 'int?', doc: "Increase or decrease to balance speed and memory usage"}
  output_basename: string
  dbnsfp: { type: 'File?', secondaryFiles: [.tbi,^.readme.txt], doc: "VEP-formatted plugin file, index, and readme file containing dbNSFP annotations" }
  dbnsfp_fields: { type: 'string?', doc: "csv string with desired fields to annotate. Use ALL to grab all" }
  merged: { type: 'boolean?', doc: "Set to true if merged cache used", default: true }
  cadd_indels: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD indel annotations" }
  cadd_snvs: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD SNV annotations" }
  run_cache_existing: { type: boolean, doc: "Run the check_existing flag for cache" }
  run_cache_af: { type: boolean, doc: "Run the allele frequency flags for cache" }
  # annotation vars
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots" }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots" }
  retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep", default: "MQ,MQ0,QSI,HotSpotAllele"}
  retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  retain_ann: { type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF" }
  add_common_fields: {type: 'boolean?', doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: true}
  bcftools_annot_columns: {type: 'string?', doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_strip_columns: {type: 'string?', doc: "csv string of columns to strip if needed to avoid conflict, i.e INFO/AF"}
  bcftools_annot_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "additional bgzipped annotation vcf file"}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration, recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP() <= 7\"), \"AF > 0.001\"]"}
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task.", default: false }
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}

  strelka2_cores: {type: 'int?', doc: "Adjust number of cores used to run strelka2", default: 16}
  extra_arg: {type: 'string?', doc: "Add special options to config file, i.e. --max-input-depth 1000"}

outputs:
  strelka2_prepass_vcf: {type: 'File', outputSource: rename_strelka_samples/reheadered_vcf}
  strelka2_protected_outputs: {type: 'File[]', outputSource: annotate/annotated_protected}
  strelka2_public_outputs: {type: 'File[]', outputSource: annotate/annotated_public}

steps:
  strelka2:
    run: ../tools/strelka2.cwl
    label: Run Strelka2
    in:
      input_tumor_aligned: input_tumor_aligned
      input_normal_aligned: input_normal_aligned
      reference: indexed_reference_fasta
      hg38_strelka_bed: hg38_strelka_bed
      manta_small_indels: manta_small_indels
      use_manta_small_indels: use_manta_small_indels
      exome_flag: exome_flag
      extra_arg: extra_arg
      cores: strelka2_cores
    out: [output_snv, output_indel]

  merge_strelka2_vcf:
    run: ../tools/gatk_mergevcfs.cwl
    label: Merge Strelka2 SNV + Indel
    in:
      input_vcfs: [strelka2/output_snv, strelka2/output_indel]
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name: tool_name
    out: [merged_vcf]

  rename_strelka_samples:
    run: ../tools/bcftools_reheader_vcf.cwl
    in:
      input_vcf: merge_strelka2_vcf/merged_vcf
      input_normal_name: input_normal_name
      input_tumor_name: input_tumor_name
    out: [reheadered_vcf]

  gatk_selectvariants_strelka2:
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select Strelka2 PASS
    in:
      input_vcf: rename_strelka_samples/reheadered_vcf
      output_basename: output_basename
      tool_name: tool_name
      mode: select_vars_mode
    out: [pass_vcf]

  annotate:
    run: ../workflow/kfdrc_annot_vcf_wf.cwl
    label: Annotate VCF
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: gatk_selectvariants_strelka2/pass_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      add_common_fields: add_common_fields
      retain_info: retain_info
      retain_fmt: retain_fmt
      retain_ann: retain_ann
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_strip_columns: bcftools_strip_columns
      bcftools_annot_vcf: bcftools_annot_vcf
      bcftools_public_filter: bcftools_public_filter
      dbnsfp: dbnsfp
      dbnsfp_fields: dbnsfp_fields
      merged: merged
      cadd_indels: cadd_indels
      cadd_snvs: cadd_snvs
      run_cache_af: run_cache_af
      run_cache_existing: run_cache_existing
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      vep_cache: vep_cache
      vep_ram: vep_ram
      vep_cores: vep_cores
      vep_buffer_size: vep_buffer_size
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      output_basename: output_basename
      tool_name: tool_name
    out: [annotated_protected, annotated_public]

$namespaces:
  sbg: https://sevenbridges.com
