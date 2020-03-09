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
  strelka2_bed: {type: File, secondaryFiles: ['.tbi']}
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
  vep_cache: {type: File, label: tar gzipped cache from ensembl/local converted cache}
  vep_ref_build: {type: ['null', string], doc: "Genome ref build used, should line up with cache.", default: "GRCh38" }
  output_basename: string
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}

outputs:
  strelka2_vep_vcf: {type: File, outputSource: vep_annot_strelka2/output_vcf}
  strelka2_vep_tbi: {type: File, outputSource: vep_annot_strelka2/output_tbi}
  strelka2_prepass_vcf: {type: File, outputSource: rename_strelka_samples/reheadered_vcf}
  strelka2_vep_maf: {type: File, outputSource: vep_annot_strelka2/output_maf}

steps:
  strelka2:
    run: ../tools/strelka2.cwl
    in:
      input_tumor_aligned: input_tumor_aligned
      input_normal_aligned: input_normal_aligned
      reference: indexed_reference_fasta
      strelka2_bed: strelka2_bed
      exome_flag: exome_flag
    out: [output_snv, output_indel]

  merge_strelka2_vcf:
    run: ../tools/gatk_mergevcfs.cwl
    label: Merge & pass filter strekla2
    in:
      input_vcfs: [strelka2/output_snv, strelka2/output_indel]
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name:
        valueFrom: ${ return "strelka2"}
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
      tool_name:
        valueFrom: ${return "strelka2"}
      mode: select_vars_mode
    out: [pass_vcf]

  vep_annot_strelka2:
    run: ../tools/vep_vcf2maf.cwl
    in:
      input_vcf: gatk_selectvariants_strelka2/pass_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name:
        valueFrom: ${return "strelka2_somatic"}
      reference: indexed_reference_fasta
      cache: vep_cache
      ref_build: vep_ref_build
    out: [output_vcf, output_tbi, output_maf, warn_txt]