cwlVersion: v1.2
class: Workflow
id: kfdrc_annot_sub_wf
doc: "This subworkflow normalizes, annotates, and adds soft filters to vcf inputs"
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  indexed_reference_fasta: {type: 'File', secondaryFiles: [.fai, ^.dict]}
  input_vcf: {type: 'File', secondaryFiles: ['.tbi'], doc: "Input vcf to annotate and soft filter"}
  input_tumor_name: string
  input_normal_name: string
  add_common_fields: {type: 'boolean', doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: false}
  # bcftools strip, if needed
  bcftools_strip_columns: {type: 'string?', doc: "csv string of columns to strip if needed to avoid conflict, i.e INFO/AF"}
  # bcftools annotate if more to do
  bcftools_annot_columns: {type: 'string?', doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF"}
  bcftools_annot_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "additional bgzipped annotation vcf file"}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration for clues"}
  # VEP-specific
  vep_ram: {type: 'int?', default: 32, doc: "In GB, may need to increase this value depending on the size/complexity of input"}
  vep_cores: {type: 'int?', default: 16, doc: "Number of cores to use. May need to increase for really large inputs"}
  vep_buffer_size: {type: 'int?', default: 5000, doc: "Increase or decrease to balance speed and memory usage"}
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache"}
  dbnsfp: { type: 'File?', secondaryFiles: [.tbi,^.readme.txt], doc: "VEP-formatted plugin file, index, and readme file containing dbNSFP annotations" }
  dbnsfp_fields: { type: 'string?', doc: "csv string with desired fields to annotate. Use ALL to grab all",
    default: 'SIFT4G_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,VEST4_score,VEST4_rankscore,MetaSVM_pred,MetaLR_pred,MetaRNN_pred,M-CAP_pred,REVEL_score,REVEL_rankscore,PrimateAI_pred,DEOGEN2_pred,BayesDel_noAF_pred,ClinPred_pred,LIST-S2_pred,Aloft_pred,fathmm-MKL_coding_pred,fathmm-XF_coding_pred,Eigen-phred_coding,Eigen-PC-phred_coding,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,TWINSUK_AC,TWINSUK_AF,ALSPAC_AC,ALSPAC_AF,UK10K_AC,UK10K_AF,gnomAD_exomes_controls_AC,gnomAD_exomes_controls_AN,gnomAD_exomes_controls_AF,gnomAD_exomes_controls_nhomalt,gnomAD_exomes_controls_POPMAX_AC,gnomAD_exomes_controls_POPMAX_AN,gnomAD_exomes_controls_POPMAX_AF,gnomAD_exomes_controls_POPMAX_nhomalt,gnomAD_genomes_flag,gnomAD_genomes_AC,gnomAD_genomes_AN,gnomAD_genomes_AF,gnomAD_genomes_nhomalt,gnomAD_genomes_POPMAX_AC,gnomAD_genomes_POPMAX_AN,gnomAD_genomes_POPMAX_AF,gnomAD_genomes_POPMAX_nhomalt,gnomAD_genomes_controls_and_biobanks_AC,gnomAD_genomes_controls_and_biobanks_AN,gnomAD_genomes_controls_and_biobanks_AF,gnomAD_genomes_controls_and_biobanks_nhomalt,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,Interpro_domain,GTEx_V8_gene,GTEx_V8_tissue'
    }
  merged: { type: 'boolean?', doc: "Set to true if merged cache used", default: true }
  run_cache_existing: { type: boolean, doc: "Run the check_existing flag for cache" }
  run_cache_af: { type: boolean, doc: "Run the allele frequency flags for cache" }
  run_stats: { type: boolean, doc: "Create stats file? Disable for speed", default: true }
  cadd_indels: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD indel annotations" }
  cadd_snvs: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD SNV annotations" }
  intervar: { type: 'File?', doc: "Intervar vcf-formatted file. See docs for custom build instructions", secondaryFiles: [.tbi] }
  # Hotspot Annotation
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task." }
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots" }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots" }
  output_basename: string
  tool_name: string
  # MAF-specific
  retain_info: { type: 'string?', doc: "csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,Hotspot`" }
  retain_fmt: { type: 'string?', doc: "csv string with FORMAT fields that you want to keep" }
  retain_ann: { type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF" }
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}

outputs:
  annotated_protected: {type: 'File[]', outputSource: rename_protected/renamed_files}
  annotated_public: {type: 'File[]', outputSource: rename_public/renamed_files}

steps:
  normalize_vcf:
    run: ../tools/normalize_vcf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: input_vcf
      output_basename: output_basename
      tool_name: tool_name
    out: [normalized_vcf]

  bcftools_strip_info:
    when: $(inputs.bcftools_strip_columns != null)
    run: ../tools/bcftools_strip_ann.cwl
    in:
      input_vcf: normalize_vcf/normalized_vcf
      output_basename: output_basename
      tool_name: tool_name
      strip_info: bcftools_strip_columns
    out: [stripped_vcf]

  add_standard_fields:
    run: ../tools/add_strelka2_fields.cwl
    when: $(inputs.run_tool_flag)
    in:
      strelka2_vcf:
        source: [bcftools_strip_info/stripped_vcf, normalize_vcf/normalized_vcf]
        pickValue: first_non_null
      run_tool_flag: add_common_fields
      tumor_name: input_tumor_name
      normal_name: input_normal_name
      output_basename: output_basename
    out: [output]

  vep_annotate_vcf:
    run: ../tools/variant_effect_predictor_105.cwl
    in:
      reference: indexed_reference_fasta
      cores: vep_cores
      ram: vep_ram
      buffer_size: vep_buffer_size
      input_vcf:
        source: [add_standard_fields/output, bcftools_strip_info/stripped_vcf, normalize_vcf/normalized_vcf]
        pickValue: first_non_null
      output_basename: output_basename
      tool_name: tool_name
      cache: vep_cache
      merged: merged
      run_cache_existing: run_cache_existing
      run_cache_af: run_cache_af
      run_stats: run_stats
      cadd_indels: cadd_indels
      cadd_snvs: cadd_snvs
      dbnsfp: dbnsfp
      dbnsfp_fields: dbnsfp_fields
      intervar: intervar
    out: [output_vcf]

  bcftools_gnomad_annotate:
    when: $(inputs.annotation_vcf != null)
    run: ../tools/bcftools_annotate.cwl
    in:
      input_vcf: vep_annotate_vcf/output_vcf
      annotation_vcf: bcftools_annot_vcf
      columns: bcftools_annot_columns
      output_basename: output_basename
      tool_name: tool_name
    out: [bcftools_annotated_vcf]

  gatk_add_soft_filter:
    run: ../tools/gatk_variant_filter.cwl
    in:
      input_vcf:
        source: [bcftools_gnomad_annotate/bcftools_annotated_vcf, vep_annotate_vcf/output_vcf]
        pickValue: first_non_null
      reference: indexed_reference_fasta
      filter_name: gatk_filter_name
      filter_expression: gatk_filter_expression
      output_basename: output_basename
      tool_name: tool_name
    out: [gatk_soft_filtered_vcf]

  hotspots_annotation:
    run: ../tools/hotspots_annotation.cwl
    in:
      input_vcf: gatk_add_soft_filter/gatk_soft_filtered_vcf
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snvs: protein_snv_hotspots
      protein_indels: protein_indel_hotspots
      output_basename: output_basename
    out: [hotspots_vcf]

  kfdrc_vcf2maf_protected:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: indexed_reference_fasta
      input_vcf: hotspots_annotation/hotspots_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name: tool_name
      retain_info: retain_info
      retain_fmt: retain_fmt
      retain_ann: retain_ann
      maf_center: maf_center
    out: [output_maf]

  hard_filter_vcf:
    run: ../tools/bcftools_filter_vcf.cwl
    in:
      input_vcf: hotspots_annotation/hotspots_vcf
      include_expression: bcftools_public_filter
      output_basename: output_basename
    out:
      [filtered_vcf]

  kfdrc_vcf2maf_public:
    run: ../tools/kf_mskcc_vcf2maf.cwl
    in:
      reference: indexed_reference_fasta
      input_vcf: hard_filter_vcf/filtered_vcf
      output_basename: output_basename
      tumor_id: input_tumor_name
      normal_id: input_normal_name
      tool_name: tool_name
      retain_info: retain_info
      retain_fmt: retain_fmt
      retain_ann: retain_ann
      maf_center: maf_center
    out: [output_maf]

  rename_protected:
    run: ../tools/generic_rename_outputs.cwl
    label: Rename Protected Outputs
    in:
      input_files:
        source: [hotspots_annotation/hotspots_vcf, kfdrc_vcf2maf_protected/output_maf]
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
        source: [hard_filter_vcf/filtered_vcf, kfdrc_vcf2maf_public/output_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: [output_basename, tool_name]
        valueFrom: "${var pub_vcf=self[0] + '.' + self[1] + '.norm.annot.public.vcf.gz'; \
        var pub_tbi=self[0] + '.' + self[1] + '.norm.annot.public.vcf.gz.tbi'; \
        var pub_maf=self[0] + '.' + self[1] + '.norm.annot.public.maf'; \
        return [pub_vcf, pub_tbi, pub_maf];}"
    out: [renamed_files]