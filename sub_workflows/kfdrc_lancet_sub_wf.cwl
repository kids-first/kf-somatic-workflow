cwlVersion: v1.0
class: Workflow
id: kfdrc_lancet_sub_wf
doc: "Lancet sub workflow, meant to be wrapped"
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  indexed_reference_fasta: {type: 'File', secondaryFiles: [.fai, ^.dict]}
  input_tumor_aligned: {type: 'File', secondaryFiles: [^.bai]}
  input_tumor_name: string
  input_normal_aligned: {type: 'File', secondaryFiles: [^.bai]}
  input_normal_name: string
  output_basename: string
  reference_dict: File
  bed_invtl_split: {type: 'File[]', doc: "Bed file intervals passed on from and outside pre-processing step"}
  ram: {type: ['null', int], default: 12, doc: "Adjust in rare circumstances in which 12 GB is not enough."}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
  window: {type: 'int', doc: "window size for lancet.  default is 600, recommend 500 for WGS, 600 for exome+"}
  padding: {type: 'int', doc: "If WGS (less likely), default 25, if exome+, recommend half window size, 0 if pure exome or panel"}
  tool_name: {type: 'string?', doc: "String to describe what tool was run as part of file name", default: "lancet_somatic"}
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache"}
  vep_ram: {type: 'int?', doc: "In GB, may need to increase this value depending on the size/complexity of input"}
  vep_cores: {type: 'int?', doc: "Number of cores to use. May need to increase for really large inputs"}
  vep_buffer_size: {type: 'int?', doc: "Increase or decrease to balance speed and memory usage"}
  dbnsfp: { type: 'File?', secondaryFiles: [.tbi,^.readme.txt], doc: "VEP-formatted plugin file, index, and readme file containing dbNSFP annotations" }
  dbnsfp_fields: { type: 'string?', doc: "csv string with desired fields to annotate. Use ALL to grab all",
    default: 'SIFT4G_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,VEST4_score,VEST4_rankscore,MetaSVM_pred,MetaLR_pred,MetaRNN_pred,M-CAP_pred,REVEL_score,REVEL_rankscore,PrimateAI_pred,DEOGEN2_pred,BayesDel_noAF_pred,ClinPred_pred,LIST-S2_pred,Aloft_pred,fathmm-MKL_coding_pred,fathmm-XF_coding_pred,Eigen-phred_coding,Eigen-PC-phred_coding,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,TWINSUK_AC,TWINSUK_AF,ALSPAC_AC,ALSPAC_AF,UK10K_AC,UK10K_AF,gnomAD_exomes_controls_AC,gnomAD_exomes_controls_AN,gnomAD_exomes_controls_AF,gnomAD_exomes_controls_nhomalt,gnomAD_exomes_controls_POPMAX_AC,gnomAD_exomes_controls_POPMAX_AN,gnomAD_exomes_controls_POPMAX_AF,gnomAD_exomes_controls_POPMAX_nhomalt,gnomAD_genomes_flag,gnomAD_genomes_AC,gnomAD_genomes_AN,gnomAD_genomes_AF,gnomAD_genomes_nhomalt,gnomAD_genomes_POPMAX_AC,gnomAD_genomes_POPMAX_AN,gnomAD_genomes_POPMAX_AF,gnomAD_genomes_POPMAX_nhomalt,gnomAD_genomes_controls_and_biobanks_AC,gnomAD_genomes_controls_and_biobanks_AN,gnomAD_genomes_controls_and_biobanks_AF,gnomAD_genomes_controls_and_biobanks_nhomalt,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,Interpro_domain,GTEx_V8_gene,GTEx_V8_tissue'
    }
  merged: { type: 'boolean?', doc: "Set to true if merged cache used", default: true }
  cadd_indels: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD indel annotations" }
  cadd_snvs: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD SNV annotations" }
  run_cache_existing: { type: boolean, doc: "Run the check_existing flag for cache" }
  run_cache_af: { type: boolean, doc: "Run the allele frequency flags for cache" }

  # annotation vars
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots" }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots" }
  retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep", default: "MS,FETS,HotSpotAllele"}
  retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  retain_ann: { type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF" }
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}
  add_common_fields: {type: 'boolean?', doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: false}
  bcftools_annot_columns: {type: 'string?', doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_strip_columns: {type: 'string?', doc: "csv string of columns to strip if needed to avoid conflict, i.e INFO/AF"}
  bcftools_annot_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "additional bgzipped annotation vcf file"}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration, recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP() <= 7\"), \"AF > 0.001\"]"}
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task.", default: false }

outputs:
  lancet_prepass_vcf: {type: 'File', outputSource: sort_merge_lancet_vcf/merged_vcf}
  lancet_protected_outputs: {type: 'File[]', outputSource: annotate/annotated_protected}
  lancet_public_outputs: {type: 'File[]', outputSource: annotate/annotated_public}

steps:
  lancet:
    run: ../tools/lancet.cwl
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge
    in:
      input_tumor_bam: input_tumor_aligned
      input_normal_bam: input_normal_aligned
      reference: indexed_reference_fasta
      bed: bed_invtl_split
      output_basename: output_basename
      window: window
      padding: padding
      ram: ram
    scatter: [bed]
    out: [lancet_vcf]

  sort_merge_lancet_vcf:
    run: ../tools/gatk_sortvcf.cwl
    label: GATK Sort & Merge lancet
    in:
      input_vcfs: lancet/lancet_vcf
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name: tool_name
    out: [merged_vcf]

  gatk_selectvariants_lancet:
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select Lancet PASS
    in:
      input_vcf: sort_merge_lancet_vcf/merged_vcf
      output_basename: output_basename
      tool_name: tool_name
      mode: select_vars_mode
    out: [pass_vcf]

  annotate:
    run: ../workflow/kfdrc_annot_vcf_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: gatk_selectvariants_lancet/pass_vcf
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
