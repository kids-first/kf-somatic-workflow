cwlVersion: v1.0
class: Workflow
id: kfdrc_strelka2_sub_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

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

outputs:
  strelka2_prepass_vcf: {type: File, outputSource: rename_strelka_samples/reheadered_vcf}
  strelka2_protected_outputs: {type: 'File[]', outputSource: rename_protected/renamed_files}
  strelka2_public_outputs: {type: 'File[]', outputSource: rename_public/renamed_files}

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
    run: ../sub_workflows/kfdrc_annot_vcf_sub_wf.cwl
    label: Annotate VCF
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: gatk_selectvariants_strelka2/pass_vcf
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
    label: Rename Protected Outputs
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
    label: Rename Public Outputs
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
