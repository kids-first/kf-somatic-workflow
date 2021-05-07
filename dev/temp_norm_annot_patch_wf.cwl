cwlVersion: v1.0
class: Workflow
id: kfdrc_norm_annot_wf
doc: "This is a temp wf to normalize and annotate a PASS vcf"
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: ['.fai', '^.dict'], sbg:suggestedValue: {
      class: File, path: 60639014357c3a53540ca7a3, name: Homo_sapiens_assembly38.fasta,
      secondaryFiles: [{class: File, path: 60639016357c3a53540ca7af, name: Homo_sapiens_assembly38.fasta},
        {class: File, path: 60639019357c3a53540ca7e7, name: Homo_sapiens_assembly38.dict}]}}
  strelka2_pass_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  mutect2_pass_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  lancet_pass_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  vardict_pass_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
  input_tumor_name: string
  input_normal_name: string
  bcftools_annot_columns: {type: string?, doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_annot_vcf: {type: File, doc: "bgzipped annotation vcf file", sbg:suggestedValue: {class: File, path: 5f50018fe4b054958bc8d2e3,
      name: af-only-gnomad.hg38.vcf.gz,  secondaryFiles: [{class: File, path: 5f50018fe4b054958bc8d2e5, name: af-only-gnomad.hg38.vcf.gz.tbi}]} }
  bcftools_public_filter: {type: string?, doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  vep_cache: {type: File, doc: "tar gzipped cache from ensembl/local converted cache",
    sbg:suggestedValue: {class: File, path: 607713829360f10e3982a425, name: homo_sapiens_vep_93_GRCh38.tar.gz}}
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots",
                      sbg:suggestedValue: [{class: File, path: 607713829360f10e3982a423, name: tert.bed}] }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots",
                      sbg:suggestedValue: [{class: File, path: 607713829360f10e3982a426, name: protein_snv_cancer_hotspots_v2.tsv}] }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots",
                      sbg:suggestedValue: [{class: File, path: 607713829360f10e3982a424, name: protein_indel_cancer_hotspots_v2.tsv}] }
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration, recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP() <= 7\"), \"AF > 0.001\"]"}
  output_basename: string
  # defaults, do not edit
  gatk_filter_name: { type: 'string[]?', doc: "GATK filter names", default: ["NORM_DP_LOW", "GNOMAD_AF_HIGH"] }
  tool_s: {type: string?, default: "strelka2_somatic"}
  tool_m: {type: string?, default: "mutect2_somatic"}
  tool_l: {type: string?, default: "lancet_somatic"}
  tool_v: {type: string?, default: "vardict_somatic"}
  add_common_fields_s: {type: boolean?, doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: true}
  add_common_fields_r: {type: boolean?, doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: false}
  strelka2_info: { type: string?, default: "MQ,MQ0,QSI,HotSpotAllele" }
  mutect2_info: { type: string?, default: "MBQ,TLOD,HotSpotAllele" }
  lancet_info: { type: string?, default: "MS,FETS,HotSpotAllele" }
  vardict_info: { type: string?, default: "MSI,MSILEN,SOR,SSF,HotSpotAllele" }

outputs:
  annotated_protected_strelka2: {type: 'File[]', outputSource: rename_strelka2_protected/renamed_files}
  annotated_public_strelka2: {type: 'File[]', outputSource: rename_strelka2_public/renamed_files}
  annotated_protected_mutect2: {type: 'File[]', outputSource: rename_mutect2_protected/renamed_files}
  annotated_public_mutect2: {type: 'File[]', outputSource: rename_mutect2_public/renamed_files}
  annotated_protected_lancet: {type: 'File[]', outputSource: rename_lancet_protected/renamed_files}
  annotated_public_lancet: {type: 'File[]', outputSource: rename_lancet_public/renamed_files}
  annotated_protected_vardict: {type: 'File[]', outputSource: rename_vardict_protected/renamed_files}
  annotated_public_vardict: {type: 'File[]', outputSource: rename_vardict_public/renamed_files}

steps:
  annotate_strelka2:
    run: ../sub_workflows/kfdrc_annot_vcf_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: strelka2_pass_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      bcftools_public_filter: bcftools_public_filter
      bcftools_annot_vcf: bcftools_annot_vcf
      bcftools_annot_columns: bcftools_annot_columns
      add_common_fields: add_common_fields_s
      retain_info: strelka2_info
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      vep_cache: vep_cache
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      output_basename: output_basename
      tool_name: tool_s
    out: [annotated_protected_vcf, annotated_protected_maf, annotated_public_vcf, annotated_public_maf]

  annotate_mutect2:
    run: ../sub_workflows/kfdrc_annot_vcf_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: mutect2_pass_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      bcftools_public_filter: bcftools_public_filter
      bcftools_annot_vcf: bcftools_annot_vcf
      bcftools_annot_columns: bcftools_annot_columns
      add_common_fields: add_common_fields_r
      retain_info: mutect2_info
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      vep_cache: vep_cache
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      output_basename: output_basename
      tool_name: tool_m
    out: [annotated_protected_vcf, annotated_protected_maf, annotated_public_vcf, annotated_public_maf]

  annotate_lancet:
    run: ../sub_workflows/kfdrc_annot_vcf_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: lancet_pass_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      bcftools_public_filter: bcftools_public_filter
      bcftools_annot_vcf: bcftools_annot_vcf
      bcftools_annot_columns: bcftools_annot_columns
      add_common_fields: add_common_fields_r
      retain_info: lancet_info
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      vep_cache: vep_cache
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      output_basename: output_basename
      tool_name: tool_l
    out: [annotated_protected_vcf, annotated_protected_maf, annotated_public_vcf, annotated_public_maf]

  annotate_vardict:
    run: ../sub_workflows/kfdrc_annot_vcf_sub_wf.cwl
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: vardict_pass_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      bcftools_public_filter: bcftools_public_filter
      bcftools_annot_vcf: bcftools_annot_vcf
      bcftools_annot_columns: bcftools_annot_columns
      add_common_fields: add_common_fields_r
      retain_info: vardict_info
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      vep_cache: vep_cache
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      output_basename: output_basename
      tool_name: tool_v
    out: [annotated_protected_vcf, annotated_protected_maf, annotated_public_vcf, annotated_public_maf]

  rename_strelka2_protected:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [annotate_strelka2/annotated_protected_vcf, annotate_strelka2/annotated_protected_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: output_basename
        valueFrom: "${var pro_vcf=self + '.strelka2_somatic.norm.annot.protected.vcf.gz'; \
        var pro_tbi=self[0] + '.strelka2_somatic.norm.annot.protected.vcf.gz.tbi'; \
        var pro_maf=self[0] + '.strelka2_somatic.norm.annot.protected.maf'; \
        return [pro_vcf, pro_tbi, pro_maf];}"
    out: [renamed_files]

  rename_strelka2_public:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [annotate_strelka2/annotated_public_vcf, annotate_strelka2/annotated_public_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: output_basename
        valueFrom: "${var pub_vcf=self + '.strelka2_somatic.norm.annot.public.vcf.gz'; \
        var pub_tbi=self[0] + '.strelka2_somatic.norm.annot.public.vcf.gz.tbi'; \
        var pub_maf=self[0] + '.strelka2_somatic.norm.annot.public.maf'; \
        return [pub_vcf, pub_tbi, pub_maf];}"
    out: [renamed_files]

  rename_mutect2_protected:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [annotate_mutect2/annotated_protected_vcf, annotate_mutect2/annotated_protected_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: output_basename
        valueFrom: "${var pro_vcf=self + '.mutect2_somatic.norm.annot.protected.vcf.gz'; \
        var pro_tbi=self[0] + '.mutect2_somatic.norm.annot.protected.vcf.gz.tbi'; \
        var pro_maf=self[0] + '.mutect2_somatic.norm.annot.protected.maf'; \
        return [pro_vcf, pro_tbi, pro_maf];}"
    out: [renamed_files]

  rename_mutect2_public:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [annotate_mutect2/annotated_public_vcf, annotate_mutect2/annotated_public_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: output_basename
        valueFrom: "${var pub_vcf=self + '.mutect2_somatic.norm.annot.public.vcf.gz'; \
        var pub_tbi=self[0] + '.mutect2_somatic.norm.annot.public.vcf.gz.tbi'; \
        var pub_maf=self[0] + '.mutect2_somatic.norm.annot.public.maf'; \
        return [pub_vcf, pub_tbi, pub_maf];}"
    out: [renamed_files]

  rename_vardict_protected:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [annotate_vardict/annotated_protected_vcf, annotate_vardict/annotated_protected_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: output_basename
        valueFrom: "${var pro_vcf=self + '.vardict_somatic.norm.annot.protected.vcf.gz'; \
        var pro_tbi=self[0] + '.vardict_somatic.norm.annot.protected.vcf.gz.tbi'; \
        var pro_maf=self[0] + '.vardict_somatic.norm.annot.protected.maf'; \
        return [pro_vcf, pro_tbi, pro_maf];}"
    out: [renamed_files]

  rename_vardict_public:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [annotate_vardict/annotated_public_vcf, annotate_vardict/annotated_public_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: output_basename
        valueFrom: "${var pub_vcf=self + '.vardict_somatic.norm.annot.public.vcf.gz'; \
        var pub_tbi=self[0] + '.vardict_somatic.norm.annot.public.vcf.gz.tbi'; \
        var pub_maf=self[0] + '.vardict_somatic.norm.annot.public.maf'; \
        return [pub_vcf, pub_tbi, pub_maf];}"
    out: [renamed_files]

  rename_lancet_protected:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [annotate_lancet/annotated_protected_vcf, annotate_lancet/annotated_protected_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: output_basename
        valueFrom: "${var pro_vcf=self + '.lancet_somatic.norm.annot.protected.vcf.gz'; \
        var pro_tbi=self[0] + '.lancet_somatic.norm.annot.protected.vcf.gz.tbi'; \
        var pro_maf=self[0] + '.lancet_somatic.norm.annot.protected.maf'; \
        return [pro_vcf, pro_tbi, pro_maf];}"
    out: [renamed_files]

  rename_lancet_public:
    run: ../tools/generic_rename_outputs.cwl
    in:
      input_files:
        source: [annotate_lancet/annotated_public_vcf, annotate_lancet/annotated_public_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: output_basename
        valueFrom: "${var pub_vcf=self + '.lancet_somatic.norm.annot.public.vcf.gz'; \
        var pub_tbi=self[0] + '.lancet_somatic.norm.annot.public.vcf.gz.tbi'; \
        var pub_maf=self[0] + '.lancet_somatic.norm.annot.public.maf'; \
        return [pub_vcf, pub_tbi, pub_maf];}"
    out: [renamed_files]


$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2
