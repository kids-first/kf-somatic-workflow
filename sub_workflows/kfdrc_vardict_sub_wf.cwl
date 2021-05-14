cwlVersion: v1.0
class: Workflow
id: kfdrc_vardict_1_7_sub_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  input_tumor_aligned: {type: File, secondaryFiles: ['^.bai']}
  input_tumor_name: string
  input_normal_aligned: {type: File, secondaryFiles: ['^.bai']}
  input_normal_name: string
  output_basename: string
  reference_dict: File
  padding: {type: ['null', int], doc: "Padding to add to input intervals, recommened 0 if intervals already padded, 150 if not", default: 150}
  bed_invtl_split: {type: 'File[]', doc: "Bed file intervals passed on from and outside pre-processing step"}
  min_vaf: {type: ['null', float], doc: "Min variant allele frequency for vardict to consider.  Recommend 0.05", default: 0.05}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
  cpus: {type: ['null', int], default: 9}
  ram: {type: ['null', int], default: 18, doc: "In GB"}
  vep_cache: {type: File, doc: "tar gzipped cache from ensembl/local converted cache"}
  vep_ref_build: {type: ['null', string], doc: "Genome ref build used, should line up with cache.", default: "GRCh38" }
  tool_name: {type: string?, doc: "String to describe what tool was run as part of file name", default: "vardict_somatic"}
  # annotation vars
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots" }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots" }
  retain_info: {type: string?, doc: "csv string with INFO fields that you want to keep", default: "MSI,MSILEN,SOR,SSF,HotSpotAllele"}
  retain_fmt: {type: string?, doc: "csv string with FORMAT fields that you want to keep"}
  add_common_fields: {type: boolean?, doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: false}
  bcftools_annot_columns: {type: string, doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_annot_vcf: {type: File, secondaryFiles: ['.tbi'], doc: "bgzipped annotation vcf file"}
  bcftools_public_filter: {type: string?, doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration, recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP() <= 7\"), \"AF > 0.001\"]"}
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task.", default: false }
  maf_center: {type: string?, doc: "Sequencing center of variant called", default: "."}

outputs:
  vardict_prepass_vcf: {type: File, outputSource: sort_merge_vardict_vcf/merged_vcf}
  vardict_protected_outputs: {type: 'File[]', outputSource: rename_protected/renamed_files}
  vardict_public_outputs: {type: 'File[]', outputSource: rename_public/renamed_files}

steps:
  vardict:
    run: ../tools/vardictjava.cwl
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge
    in:
      input_tumor_bam: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_bam: input_normal_aligned
      input_normal_name: input_normal_name
      padding: padding
      min_vaf: min_vaf
      cpus: cpus
      ram: ram
      reference: indexed_reference_fasta
      bed: bed_invtl_split
      output_basename: output_basename
    scatter: [bed]
    out: [vardict_vcf]

  sort_merge_vardict_vcf:
    run: ../tools/gatk_sortvcf.cwl
    label: GATK Sort & merge vardict
    in:
      input_vcfs: vardict/vardict_vcf
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name: tool_name
    out: [merged_vcf]

  bcbio_filter_fp_somatic:
    run: ../tools/bcbio_filter_vardict_somatic.cwl
    in:
      input_vcf: sort_merge_vardict_vcf/merged_vcf
      output_basename: output_basename
    out: [filtered_vcf] 

  gatk_selectvariants_vardict:
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select Vardict PASS
    in:
      input_vcf: bcbio_filter_fp_somatic/filtered_vcf
      output_basename: output_basename
      tool_name: tool_name
      mode: select_vars_mode
    out: [pass_vcf]

  annotate:
    run: ../sub_workflows/kfdrc_annot_vcf_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: gatk_selectvariants_vardict/pass_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      add_common_fields: add_common_fields
      retain_info: retain_info
      retain_fmt: retain_fmt
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_annot_vcf: bcftools_annot_vcf
      bcftools_public_filter: bcftools_public_filter
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      output_basename: output_basename
      tool_name: tool_name
    out: [annotated_protected_vcf, annotated_protected_maf, annotated_public_vcf, annotated_public_maf]

  rename_protected:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [annotate/annotated_protected_vcf, annotate/annotated_protected_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: [output_basename, tool_name]
        valueFrom: "${var pro_vcf=self[0] + '.' + self[1] + '.norm.annot.protected.vcf.gz'; \
        var pro_tbi=self[0] + '.' + self[1] + '.norm.annot.protected.vcf.gz.tbi'; \
        var pro_maf=self[0] + '.' + self[1] + '.norm.annot.protected.maf'; \
        return [pro_vcf, pro_tbi, pro_maf];}"
    out: [renamed_files]

  rename_public:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [annotate/annotated_public_vcf, annotate/annotated_public_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: [output_basename, tool_name]
        valueFrom: "${var pub_vcf=self[0] + '.' + self[1] + '.norm.annot.public.vcf.gz'; \
        var pub_tbi=self[0] + '.' + self[1] + '.norm.annot.public.vcf.gz.tbi'; \
        var pub_maf=self[0] + '.' + self[1] + '.norm.annot.public.maf'; \
        return [pub_vcf, pub_tbi, pub_maf];}"
    out: [renamed_files]

$namespaces:
  sbg: https://sevenbridges.com
