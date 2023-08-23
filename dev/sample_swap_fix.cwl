cwlVersion: v1.0
class: Workflow
id: sample_swap_fix 
doc: |
  Fix sample swap error.

  Input:
  - Mutect2, Lancet, and Vardict prepass VCFs with swapped samples
  - Strelka2 Annotated Protected VCF
requirements:
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement

inputs:
  indexed_reference_fasta: {type: 'File', secondaryFiles: ['.fai', '^.dict'], "sbg:suggestedValue": {
      class: File, path: 60639014357c3a53540ca7a3, name: Homo_sapiens_assembly38.fasta,
      secondaryFiles: [{class: File, path: 60639016357c3a53540ca7af, name: Homo_sapiens_assembly38.fasta},
        {class: File, path: 60639019357c3a53540ca7e7, name: Homo_sapiens_assembly38.dict}]}}
  strelka2_protected_vcf: {type: 'File', secondaryFiles: ['.tbi']}
  bad_prepass_mutect2_vcf: {type: 'File', secondaryFiles: ['.tbi']}
  bad_prepass_lancet_vcf: {type: 'File', secondaryFiles: ['.tbi']}
  bad_prepass_vardict_vcf: {type: 'File', secondaryFiles: ['.tbi']}
  cram: {type: 'File', secondaryFiles: ['.crai'], doc: "Tumor cram recommended for\
      \ MQ score calculation"}
  input_tumor_name: string
  input_normal_name: string
  output_basename_somatic: string
  output_basename_consensus: string
  ncallers: {type: 'int?', doc: "Optional number of callers required for consensus\
      \ [2]", default: 2}
  hotspot_source: {type: 'string?', doc: "Optional description of hotspot definition\
      \ source"}
  contig_bed: {type: 'File?', doc: "Optional BED file containing names of target contigs\
      \ / chromosomes"}
  consensus_ram: {type: 'int?', doc: "Set min memory in GB for consensus merge step",
    default: 3}
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache",
    "sbg:suggestedValue": {class: File, path: 6332f8e47535110eb79c794f, name: homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz}}
  dbnsfp: {type: 'File?', secondaryFiles: [.tbi, ^.readme.txt], doc: "VEP-formatted\
      \ plugin file, index, and readme file containing dbNSFP annotations"}
  dbnsfp_fields: {type: 'string?', doc: "csv string with desired fields to annotate\
      \ if dbnsfp provided. Use ALL to grab all"}
  merged: {type: 'boolean?', doc: "Set to true if merged cache used", default: true}
  cadd_indels: {type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin\
      \ file and index containing CADD indel annotations"}
  cadd_snvs: {type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file\
      \ and index containing CADD SNV annotations"}
  run_cache_existing: {type: 'boolean?', doc: "Run the check_existing flag for cache"}
  run_cache_af: {type: 'boolean?', doc: "Run the allele frequency flags for cache"}

  # MAF-specific params
  lancet_retain_info: {type: 'string?', doc: "csv string with INFO fields that you\
      \ want to keep, i.e. for lancet `MS,FETS,HotSpotAllele`", default: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,MBQ,TLOD,HotSpotAllele"}
  lancet_retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you\
      \ want to keep"}
  lancet_retain_ann: {type: 'string?', doc: "csv string of annotations (within the\
      \ VEP CSQ/ANN) to retain as extra columns in MAF", default: "HGVSg"}
  mutect2_retain_info: {type: 'string?', doc: "csv string with INFO fields that you\
      \ want to keep, i.e. for mutect2 `MBQ,TLOD,HotSpotAllele`", default: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,MBQ,TLOD,HotSpotAllele"}
  mutect2_retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you\
      \ want to keep"}
  mutect2_retain_ann: {type: 'string?', doc: "csv string of annotations (within the\
      \ VEP CSQ/ANN) to retain as extra columns in MAF", default: "HGVSg"}
  vardict_retain_info: {type: 'string?', doc: "csv string with INFO fields that you\
      \ want to keep, i.e. for consensus `MSI,MSILEN,SOR,SSF,HotSpotAllele`", default: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,MSI,MSILEN,SOR,SSF,HotSpotAllele"}
  vardict_retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you\
      \ want to keep"}
  vardict_retain_ann: {type: 'string?', doc: "csv string of annotations (within the\
      \ VEP CSQ/ANN) to retain as extra columns in MAF", default: "HGVSg"}
  consensus_retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to\
      \ keep", default: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,MQ,MQ0,CAL,HotSpotAllele"}
  consensus_retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want\
      \ to keep"}
  consensus_retain_ann: {type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN)\
      \ to retain as extra columns in MAF", default: "HGVSg"}

  # annotation vars
  genomic_hotspots: {type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing\
      \ hg38 genomic positions corresponding to hotspots", "sbg:suggestedValue": [
      {class: File, path: 607713829360f10e3982a423, name: tert.bed}]}
  protein_snv_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited\
      \ file(s) containing protein names and amino acid positions corresponding to\
      \ hotspots", "sbg:suggestedValue": [{class: File, path: 645919782fe81458768c552c,
        name: protein_snv_cancer_hotspots_v2.ENS105_liftover.tsv}]}
  protein_indel_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited\
      \ file(s) containing protein names and amino acid position ranges corresponding\
      \ to hotspots", "sbg:suggestedValue": [{class: File, path: 645919782fe81458768c552d,
        name: protein_indel_cancer_hotspots_v2.ENS105_liftover.tsv}]}
  bcftools_annot_columns: {type: 'string?', doc: "csv string of columns from annotation\
      \ to port into the input vcf, i.e INFO/AF", default: "INFO/gnomad_3_1_1_AC:=INFO/AC,INFO/gnomad_3_1_1_AN:=INFO/AN,INFO/gnomad_3_1_1_AF:=INFO/AF,INFO/gnomad_3_1_1_nhomalt:=INFO/nhomalt,INFO/gnomad_3_1_1_AC_popmax:=INFO/AC_popmax,INFO/gnomad_3_1_1_AN_popmax:=INFO/AN_popmax,INFO/gnomad_3_1_1_AF_popmax:=INFO/AF_popmax,INFO/gnomad_3_1_1_nhomalt_popmax:=INFO/nhomalt_popmax,INFO/gnomad_3_1_1_AC_controls_and_biobanks:=INFO/AC_controls_and_biobanks,INFO/gnomad_3_1_1_AN_controls_and_biobanks:=INFO/AN_controls_and_biobanks,INFO/gnomad_3_1_1_AF_controls_and_biobanks:=INFO/AF_controls_and_biobanks,INFO/gnomad_3_1_1_AF_non_cancer:=INFO/AF_non_cancer,INFO/gnomad_3_1_1_primate_ai_score:=INFO/primate_ai_score,INFO/gnomad_3_1_1_splice_ai_consequence:=INFO/splice_ai_consequence"}
  bcftools_annot_vcf: {type: 'File', doc: "bgzipped annotation vcf file", "sbg:suggestedValue": {
      class: File, path: 6324ef5ad01163633daa00d8, name: gnomad_3.1.1.vwb_subset.vcf.gz, secondaryFiles: [{class: File, path: 6324ef5ad01163633daa00d7, name: gnomad_3.1.1.vwb_subset.vcf.gz.tbi}]}}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to\
      \ create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to\
      \ add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to\
      \ establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration,\
      \ recommend: `vc.getGenotype('inputs.input_normal_name').getDP()\
      \ <= 7)`, `gnomad_3_1_1_AF > 0.001`]"}
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}
  custom_enst: { type: 'File?', doc: "Use a file with ens tx IDs for each gene to override VEP PICK", "sbg:suggestedValue": {class: File, path: 6480c8a61dfc710d24a3a368,
        name: kf_isoform_override.tsv} }

  # annotation vars
  bcftools_strip_columns: {type: 'string?', doc: "csv string of columns to strip if\
      \ needed to avoid conflict, i.e INFO/AF"}

  # Theta2
  cnvkit_calls: { type: 'File', doc: "cnvkit_calls from the somatic workflow" }
  cnvkit_cnn_output: { type: 'File', doc: "cnvkit_cnn_output from the somatic workflow" }
  combined_include_expression: {type: 'string?', doc: "Theta2 Purity value: Filter\
      \ expression if vcf has non-PASS combined calls, use as-needed, default set\
      \ for VarDict Java for VarDict", default: FILTER="PASS" && (INFO/STATUS="Germline"
      | INFO/STATUS="StrongSomatic")}
  combined_exclude_expression: {type: 'string?', doc: "Theta2 Purity value: Filter\
      \ expression if vcf has non-PASS combined calls, use as-needed"}
  min_theta2_frac: {type: 'float?', default: 0.01, doc: "Minimum fraction of genome\
      \ with copy umber alterations.  Default is 0.05, recommend 0.01"}

outputs:
  lancet_prepass_vcf: {type: 'File', outputSource: rename_lancet_vcf_samples/reheadered_vcf}
  mutect2_prepass_vcf: {type: 'File', outputSource: rename_mutect2_vcf_samples/reheadered_vcf}
  vardict_prepass_vcf: {type: 'File', outputSource: rename_vardict_vcf_samples/reheadered_vcf}
  theta2_calls: {type: 'File?', outputSource: run_theta2_purity/theta2_adjusted_cns}
  theta2_seg: {type: 'File?', outputSource: run_theta2_purity/theta2_adjusted_seg}
  theta2_subclonal_results: {type: ['null', 'File[]'], outputSource: expression_flatten_subclonal_results/output}
  theta2_subclonal_cns: {type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclonal_cns}
  theta2_subclone_seg: {type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclone_seg}
  lancet_protected_outputs: {type: 'File[]', outputSource: annotate_lancet/annotated_protected}
  lancet_public_outputs: {type: 'File[]', outputSource: annotate_lancet/annotated_public}
  mutect2_protected_outputs: {type: 'File[]', outputSource: annotate_mutect2/annotated_protected}
  mutect2_public_outputs: {type: 'File[]', outputSource: annotate_mutect2/annotated_public}
  vardict_protected_outputs: {type: 'File[]', outputSource: annotate_vardict/annotated_protected}
  vardict_public_outputs: {type: 'File[]', outputSource: annotate_vardict/annotated_public}
  consensus_protected_outputs: {type: 'File[]', outputSource: annotate_consensus/annotated_protected}
  consensus_public_outputs: {type: 'File[]', outputSource: annotate_consensus/annotated_public}

steps:
  rename_lancet_vcf_samples:
    run: ../tools/bcftools_reheader_samples_index.cwl
    in:
      input_vcf: bad_prepass_lancet_vcf
      output_filename:
        source: output_basename_somatic
        valueFrom: |
          $(self).lancet_somatic.merged.reheadered.vcf.gz
      new_normal_name: input_normal_name
      new_tumor_name: input_tumor_name
      old_normal_name: input_tumor_name
      old_tumor_name: input_normal_name
      tbi:
        valueFrom: |
          $(1 == 1)
    out: [reheadered_vcf]
  rename_mutect2_vcf_samples:
    run: ../tools/bcftools_reheader_samples_index.cwl
    in:
      input_vcf: bad_prepass_mutect2_vcf
      output_filename:
        source: output_basename_somatic
        valueFrom: |
          $(self).mutect2_filtered.merged.reheadered.vcf.gz
      new_normal_name: input_normal_name
      new_tumor_name: input_tumor_name
      old_normal_name: input_tumor_name
      old_tumor_name: input_normal_name
      tbi:
        valueFrom: |
          $(1 == 1)
    out: [reheadered_vcf]
  rename_vardict_vcf_samples:
    run: ../tools/bcftools_reheader_samples_index.cwl
    in:
      input_vcf: bad_prepass_vardict_vcf 
      output_filename:
        source: output_basename_somatic
        valueFrom: |
          $(self).vardict_somatic.merged.reheadered.vcf.gz
      new_normal_name: input_normal_name
      new_tumor_name: input_tumor_name
      old_normal_name: input_tumor_name
      old_tumor_name: input_normal_name
      tbi:
        valueFrom: |
          $(1 == 1)
    out: [reheadered_vcf]
  gatk_selectvariants_lancet:
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select Lancet PASS
    in:
      input_vcf: rename_lancet_vcf_samples/reheadered_vcf 
      output_basename: output_basename_somatic
      tool_name:
        valueFrom: "lancet_somatic"
      mode:
        valueFrom: "gatk"
    out: [pass_vcf]
  gatk_selectvariants_mutect2:
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select PASS
    in:
      input_vcf: rename_mutect2_vcf_samples/reheadered_vcf
      output_basename: output_basename_somatic
      tool_name:
        valueFrom: "mutect2_somatic"
      mode:
        valueFrom: "gatk"
    out: [pass_vcf]
  bcbio_filter_fp_somatic:
    run: ../tools/bcbio_filter_vardict_somatic.cwl
    in:
      input_vcf: rename_vardict_vcf_samples/reheadered_vcf
      output_basename: output_basename_somatic
    out: [filtered_vcf]
  gatk_selectvariants_vardict:
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select Vardict PASS
    in:
      input_vcf: bcbio_filter_fp_somatic/filtered_vcf
      output_basename: output_basename_somatic
      tool_name:
        valueFrom: "vardict_somatic"
      mode:
        valueFrom: "gatk"
    out: [pass_vcf]
  annotate_lancet:
    run: ../workflow/kfdrc_annot_vcf_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: gatk_selectvariants_lancet/pass_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      add_common_fields:
        valueFrom: $(1 == 0)
      retain_info: lancet_retain_info
      retain_fmt: lancet_retain_fmt
      retain_ann: lancet_retain_ann
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
      disable_hotspot_annotation:
        valueFrom: $(1 == 0)
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      custom_enst: custom_enst
      output_basename: output_basename_somatic
      tool_name:
        valueFrom: "lancet_somatic"
    out: [annotated_protected, annotated_public]
  annotate_mutect2:
    run: ../workflow/kfdrc_annot_vcf_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: gatk_selectvariants_mutect2/pass_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      add_common_fields:
        valueFrom: $(1 == 0)
      retain_info: mutect2_retain_info
      retain_fmt: mutect2_retain_fmt
      retain_ann: mutect2_retain_ann
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
      disable_hotspot_annotation:
        valueFrom: $(1 == 0)
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      custom_enst: custom_enst
      output_basename: output_basename_somatic
      tool_name:
        valueFrom: "mutect2_somatic"
    out: [annotated_protected, annotated_public]
  annotate_vardict:
    run: ../workflow/kfdrc_annot_vcf_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: gatk_selectvariants_vardict/pass_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      add_common_fields:
        valueFrom: $(1 == 0)
      retain_info: vardict_retain_info
      retain_fmt: vardict_retain_fmt
      retain_ann: vardict_retain_ann
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
      disable_hotspot_annotation:
        valueFrom: $(1 == 0)
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      custom_enst: custom_enst
      output_basename: output_basename_somatic
      tool_name:
        valueFrom: "vardict_somatic"
    out: [annotated_protected, annotated_public]

  # Theta2

  run_theta2_purity:
    run: ../sub_workflows/kfdrc_run_theta2_sub_wf.cwl
    in:
      tumor_cns: cnvkit_calls
      reference_cnn: cnvkit_cnn_output
      tumor_sample_name: input_tumor_name
      normal_sample_name: input_normal_name
      paired_vcf: rename_vardict_vcf_samples/reheadered_vcf
      combined_include_expression: combined_include_expression
      combined_exclude_expression: combined_exclude_expression
      min_theta2_frac: min_theta2_frac
      output_basename: output_basename_somatic
    out: [theta2_adjusted_cns, theta2_adjusted_seg, theta2_subclonal_results, theta2_subclonal_cns,
      theta2_subclone_seg]

  expression_flatten_subclonal_results:
    run: ../tools/expression_flatten_file_list.cwl
    in:
      input_list: run_theta2_purity/theta2_subclonal_results
    out: [output]

  # CONSENSUS

  prep_mnp_variants:
    run: ../tools/prep_mnp_variants.cwl
    in:
      strelka2_vcf: strelka2_protected_vcf
      other_vcfs:
        source: [annotate_lancet/annotated_protected, annotate_mutect2/annotated_protected, annotate_vardict/annotated_protected]
        valueFrom: |
          ${
            var lancet_vcf = self[0].filter(function(e) { return e.nameext == ".gz" })[0];
            lancet_vcf.secondaryFiles = self[0].filter(function(e) { return e.nameext == ".tbi" });
            var mutect2_vcf = self[1].filter(function(e) { return e.nameext == ".gz" })[0];
            mutect2_vcf.secondaryFiles = self[1].filter(function(e) { return e.nameext == ".tbi" });
            var vardict_vcf = self[2].filter(function(e) { return e.nameext == ".gz" })[0];
            vardict_vcf.secondaryFiles = self[2].filter(function(e) { return e.nameext == ".tbi" });
            return [lancet_vcf, mutect2_vcf, vardict_vcf];
          }
      output_basename: output_basename_consensus
    out: [output_vcfs]

  consensus_merge:
    run: ../tools/consensus_merge.cwl
    in:
      strelka2_vcf:
        source: prep_mnp_variants/output_vcfs
        valueFrom: '$(self[0])'
      mutect2_vcf:
        source: annotate_mutect2/annotated_protected
        valueFrom: |
          ${
            var vcf = self.filter(function(e) { return e.nameext == ".gz" })[0];
            vcf.secondaryFiles = self.filter(function(e) { return e.nameext == ".tbi" });
            return vcf;
          }
      lancet_vcf:
        source: annotate_lancet/annotated_protected
        valueFrom: |
          ${
            var vcf = self.filter(function(e) { return e.nameext == ".gz" })[0];
            vcf.secondaryFiles = self.filter(function(e) { return e.nameext == ".tbi" });
            return vcf;
          }
      vardict_vcf:
        source: annotate_vardict/annotated_protected
        valueFrom: |
          ${
            var vcf = self.filter(function(e) { return e.nameext == ".gz" })[0];
            vcf.secondaryFiles = self.filter(function(e) { return e.nameext == ".tbi" });
            return vcf;
          }
      cram: cram
      ncallers: ncallers
      ram: consensus_ram
      reference: indexed_reference_fasta
      output_basename: output_basename_consensus
      hotspot_source: hotspot_source
      contig_bed: contig_bed
    out: [output]

  annotate_consensus:
    run: ../workflow/kfdrc_annot_vcf_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: consensus_merge/output
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      add_common_fields:
        valueFrom: $(1 == 0)
      retain_info: consensus_retain_info
      retain_fmt: consensus_retain_fmt
      retain_ann: consensus_retain_ann
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
      disable_hotspot_annotation:
        valueFrom: $(1 == 1)
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      custom_enst: custom_enst
      output_basename: output_basename_consensus
      tool_name:
        valueFrom: "consensus_somatic"
    out: [annotated_protected, annotated_public]

$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 3
"sbg:license": Apache License 2.0
"sbg:publisher": KFDRC
"sbg:links":
- id: 'https://github.com/kids-first/kf-somatic-workflow/releases/tag/v4.3.5'
  label: github-release
