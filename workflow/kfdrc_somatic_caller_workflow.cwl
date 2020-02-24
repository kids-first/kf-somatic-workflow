cwlVersion: v1.0
class: Workflow
doc: >-
  This is the Kids First Data Resource Center (DRC) Whole Genome Sequencing (WGS) Somatic Workflow, which includes somatic variant calling, copy number variation (CNV), and structural variant (SV) calls.
  
  
  ![data service logo](https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcS9BnbvIsTkK3QlSGMDvlgu0tZQJ1q4crMvA-S3fcWfIq6y2d2Y)
  
  
  This workflow takes aligned cram input and performs somatic variant calling using Strelka2 and Mutect2, CNV estimation using ControlFreeC, and SV calls using Manta.
  Somatic variant and SV call results are annoated using Variant Effect Predictor, with the Memorial Sloane Kettering Cancer Center (MSKCC) vcf2maf wrapper.
  
  ### Recent updates
  
  As of February 24, 2020, this workflow has been updated to make b allele (germline call) input file for copy number truly optional.
  Also, some [GATK-recommended](https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set) filters are applied to input file, plus a min DP 10 requirement, when given
  A brief description of what this file is and a way to generate it is found in the Tips to Run section.
  Also, vcf2maf version has been updated as the previous version had bug handling Strelka2 input.

  ### Somatic Variant Calling:
  
  [Strelka2](https://github.com/Illumina/strelka) v2.9.3 calls single nucelotide variants (SNVS) and insertions/deletions (INDELS).
  [Mutect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php) v4.1.10 from the Broad institute also calls SNVS and INDELS.
  Each caller has a different approach to variant calling, and together one can glean confident results.  Strelka2 is run with default settings, similarly Mutect2 following Broad Best Practices, as of this [workflow](https://github.com/broadinstitute/gatk/blob/4.1.1.0/scripts/mutect2_wdl/mutect2.wdl).
  Futhermore, each tool's results, in variant call format (vcf), are filtered on the `PASS` flag. 
  The pre-`PASS` filtered results can still be obtained from the workflow in the eent the user wishes to keep some calls that failed `PASS` criteria.

  ### CNV Estimation:
  
  [ControlFreeC](https://github.com/BoevaLab/FREEC) v11.6 is used for CNV estimation.
  The tool portion of the workflow is a port from the Seven Bridges Genomics team, with a slight tweak in image outputs.
  Also, the workflow wrapper limits what inputs and outputs are used based on our judgement of utility.
  Outputs include raw ratio calls, copy number cals with p values assigned, b allele frequency data, as well as copy number and b allele frequency plots.

  ### SV Calls:
  
  [Manta](https://github.com/Illumina/manta) v1.4.0 is used to call SVs. Output is also in vcf format, with calls filtered on `PASS`.
  Default settings are used at run time.

  ### Variant Annotation
  
  [Variant Effect Predictor](https://useast.ensembl.org/info/docs/tools/vep/index.html) release 93, wrapped by [vcf2maf](https://github.com/mskcc/vcf2maf) v1.6.16 is used to annotate somatic variant and SV calls.
  Both the annotated vcf and maf file are made available.

  ### Tips To Run:

  1) For input cram files, be sure to have indexed them beforehand as well.
  
  2) For ControlFreeC, it is highly recommended that you supply a vcf file with germline calls (`b_allele` file input), GATK Haplotype caller (can use [this workflow](https://github.com/kids-first/kf-jointgenotyping-workflow/blob/master/workflow/kfdrc_single_sample_genotype_basic.cwl) from our git repo) recommended.
  This input is optional, but has been shown to increase copy number call accuracy.
  Please also make sure the index for this file is available.
  Also, a range of input ploidy possibilities for the inputs are needed.  You can simply use `2`, or put in a range, as an array, like 2, 3, 4.
  
  3) As a cavatica app, default references for hg38 are already pre-populated, as well as some default settings - i.e., number of threads, coefficient of variation input for ControlFreec, and `PASS` filter tool mode.

  4) What is `select_vars_mode` you ask?  On occasion, using GATK's `SelectVariants` tool will fail, so a simple `grep` mode on `PASS` can be used instead.
  Related, `bcftools_filter_vcf` is built in as a convenience in case youe b allele frequency file has not been filtered on `PASS`.
  You can use the `include_expression` `Filter="PASS"` to achieve this.
  
  5) Suggested inputs are:

    - reference_fasta: Homo_sapiens_assembly38.fasta # Genome reference fasta file
  
    - reference_dict: Homo_sapiens_assembly38.dict # Sequence dictionary created using GATK CreateSequenceDictionary from genome fasta file
  
    - wgs_calling_interval_list: wgs_canonical_calling_regions.hg38.interval_lis # GATK interval calling list, with chr 1-22, X,Y,M
  
    - af_only_gnomad_vcf: af-only-gnomad.hg38.vcf.gz # Broad GATK gnomad reference file
  
    - exac_common_vcf: small_exac_common_3.hg38.vcf.gz # Broad GATK exac reference file
  
    - hg38_strelka_bed: hg38_strelka.bed.gz # bgzipped chromosome bed file, chr 1-22, X, Y ,M
  
    - threads: 16
  
    - chr_len: hs38_chr.len # file with chromosome lengths
  
    - coeff_var: 0.05
  
    - contamination_adjustment: FALSE
  
    - vep_cache: homo_sapiens_vep_93_GRCh38_convert_cache.tar.gz # tar gzipped cache from ensembl/local converted cache
  
   - include_expression: `Filter="PASS"`

  ### Links/Resources:
  
  The related Github branch for this app is located [here](https://github.com/kids-first/kf-somatic-workflow/tree/mb-publish).
  From here you can pull the project, and modify tool params, inputs, outputs, etc. if you are a developer or a an ambitious researcher.

sbg:license: "Apache License 2.0"
sbg:categories: ["SOMATIC"]
id: kfdrc_somatic_wf
label: Kids First DRC Somatic WGS Workflow
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  reference_fasta: {type: File, doc: "Genome reference fasta file", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe6c', name: 'Homo_sapiens_assembly38.fasta'}}
  reference_dict: {type: File, doc: "Sequence dictionary created using GATK CreateSequenceDictionary from genome fasta file", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe71', name: 'Homo_sapiens_assembly38.dict'}}
  wgs_calling_interval_list: {type: File, doc: "GATK interval calling list", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe6e', name: 'wgs_canonical_calling_regions.hg38.interval_list'}}
  af_only_gnomad_vcf: {type: File, doc: "Broad GATK gnomad reference file", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe70', name: 'af-only-gnomad.hg38.vcf.gz'}}
  exac_common_vcf: {type: File, doc: "Broad GATK exac reference file", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe72', name: 'small_exac_common_3.hg38.vcf.gz'}}
  hg38_strelka_bed: {type: File, doc: "bgzipped chromosome bed file", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe6d', name: 'hg38_strelka.bed.gz'}}
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
  cfree_threads: {type: ['null', int], doc: "For ControlFreeC.  Recommend 16 max, as I/O gets saturated after that losing any advantage.", default: 16}
  # exome_flag: {type: ['null', string], doc: "insert 'Y' if exome mode"}
  # capture_regions: {type: ['null', File], doc: "If not WGS, provide this bed file"}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
  vep_cache: {type: File, doc: "tar gzipped cache from ensembl/local converted cache", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe6f', name: 'homo_sapiens_vep_93_GRCh38_convert_cache.tar.gz'}}
  output_basename: string
  cfree_ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try"}
  cfree_mate_orientation_sample: {type: ['null', {type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  cfree_mate_orientation_control: {type: ['null', {type: enum, name: mate_orientation_control, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  b_allele: {type: ['null', File], doc: "germline calls, needed for BAF, optional.  Calls from GATK haplotype caller are recommended.  Tool will prefilter germline PASS and apply GATK-recommended filter options"}
  cfree_chr_len: {type: File, doc: "TSV with chromsome names and lengths. Limit to chromosomes you actually want analyzed in ControlFreeC", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe73', name: 'hs38_chr.len'}}
  cfree_coeff_var: {type: ['null', float], default: 0.05, doc: "Coefficient of variation to set window size.  Default 0.05 recommended"}
  cfree_contamination_adjustment: {type: ['null', boolean], doc: "TRUE or FALSE to have ControlFreec estimate normal contam"}
  cfree_sex: {type: ['null', {type: enum, name: sex, symbols: ["XX", "XY"] }], doc: "If known, XX for female, XY for male"}

outputs:
  strelka2_vep_vcf: {type: File, outputSource: run_strelka2/strelka2_vep_vcf}
  strelka2_vep_tbi: {type: File, outputSource: run_strelka2/strelka2_vep_tbi}
  strelka2_prepass_vcf: {type: File, outputSource: run_strelka2/strelka2_prepass_vcf}
  strelka2_vep_maf: {type: File, outputSource: run_strelka2/strelka2_vep_maf}
  mutect2_vep_vcf: {type: File, outputSource: run_mutect2/mutect2_vep_vcf}
  mutect2_vep_tbi: {type: File, outputSource: run_mutect2/mutect2_vep_tbi}
  mutect2_prepass_vcf: {type: File, outputSource: run_mutect2/mutect2_filtered_vcf}
  mutect2_vep_maf: {type: File, outputSource: run_mutect2/mutect2_vep_maf}
  manta_pass_vcf: {type: File, outputSource: run_manta/manta_pass_vcf}
  manta_prepass_vcf: {type: File, outputSource: run_manta/manta_prepass_vcf}
  ctrlfreec_pval: {type: File, outputSource: run_controlfreec/ctrlfreec_pval}
  ctrlfreec_config: {type: File, outputSource: run_controlfreec/ctrlfreec_config}
  ctrlfreec_pngs: {type: 'File[]', outputSource: run_controlfreec/ctrlfreec_pngs}
  ctrlfreec_bam_ratio: {type: File, outputSource: run_controlfreec/ctrlfreec_bam_ratio}
  ctrlfreec_bam_seg: {type: File, outputSource: run_controlfreec/ctrlfreec_bam_seg}
  ctrlfreec_baf: {type: File?, outputSource: run_controlfreec/ctrlfreec_baf}
  ctrlfreec_info: {type: File, outputSource: run_controlfreec/ctrlfreec_info}

steps:

  index_references:
    run: ../tools/index_references.cwl
    in:
      input_fasta_file: reference_fasta
      af_only_gnomad_vcf: af_only_gnomad_vcf
      exac_common_vcf: exac_common_vcf
      hg38_strelka_bed: hg38_strelka_bed
    out:
      [indexed_reference_fasta, indexed_af_only_gnomad_vcf, indexed_exac_common_vcf, indexed_hg38_strelka_bed, reference_fai]

  gatk_filter_germline:
    run: ../tools/gatk_filter_germline_variant.cwl
    in:
      input_vcf: b_allele
      reference_fasta: index_references/indexed_reference_fasta
      output_basename: output_basename
    out:
      [filtered_vcf, filtered_pass_vcf]

  samtools_tumor_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_tumor_aligned
      threads:
        valueFrom: ${return 36}
      reference: index_references/indexed_reference_fasta
    out: [bam_file]

  samtools_normal_cram2bam:
    run: ../tools/samtools_cram2bam.cwl
    in:
      input_reads: input_normal_aligned
      threads:
        valueFrom: ${return 36}
      reference: index_references/indexed_reference_fasta
    out: [bam_file]

  run_controlfreec:
    run: ../sub_workflows/kfdrc_controlfreec_sub_wf.cwl
    in:
      input_tumor_aligned: samtools_tumor_cram2bam/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_normal_cram2bam/bam_file
      threads: cfree_threads
      output_basename: output_basename
      ploidy: cfree_ploidy
      mate_orientation_sample: cfree_mate_orientation_sample
      mate_orientation_control: cfree_mate_orientation_control
      indexed_reference_fasta: index_references/indexed_reference_fasta
      reference_fai: index_references/reference_fai
      b_allele: gatk_filter_germline/filtered_pass_vcf
      chr_len: cfree_chr_len
      coeff_var: cfree_coeff_var
      contamination_adjustment: cfree_contamination_adjustment
      cfree_sex: cfree_sex
    out:
      [ctrlfreec_pval, ctrlfreec_config, ctrlfreec_pngs, ctrlfreec_bam_ratio, ctrlfreec_bam_seg, ctrlfreec_baf, ctrlfreec_info]

  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: wgs_calling_interval_list
      reference_dict: reference_dict
      # exome_flag: exome_flag
      scatter_ct:
        valueFrom: ${return 50}
      bands:
        valueFrom: ${return 80000000}
    out: [output]

  run_strelka2:
    run: ../sub_workflows/kfdrc_strelka2_sub_wf.cwl
    in:
      indexed_reference_fasta: index_references/indexed_reference_fasta
      reference_dict: reference_dict
      hg38_strelka_bed: index_references/indexed_hg38_strelka_bed
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      exome_flag:
        valueFrom: ${return "N";}
      vep_cache: vep_cache
      output_basename: output_basename
      select_vars_mode: select_vars_mode
    out:
      [strelka2_vep_vcf, strelka2_vep_tbi, strelka2_prepass_vcf, strelka2_vep_maf]

  run_manta:
    run: ../sub_workflows/kfdrc_manta_sub_wf.cwl
    in:
      indexed_reference_fasta: index_references/indexed_reference_fasta
      reference_dict: reference_dict
      hg38_strelka_bed: index_references/indexed_hg38_strelka_bed
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      vep_cache: vep_cache
      output_basename: output_basename
      select_vars_mode: select_vars_mode
    out:
      [manta_prepass_vcf, manta_pass_vcf]


  run_mutect2:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge
    run: ../sub_workflows/kfdrc_mutect2_sub_wf.cwl
    in:
      indexed_reference_fasta: index_references/indexed_reference_fasta
      reference_dict: reference_dict
      bed_invtl_split: gatk_intervallisttools/output
      af_only_gnomad_vcf: index_references/indexed_af_only_gnomad_vcf
      exac_common_vcf: index_references/indexed_exac_common_vcf
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      exome_flag:
        valueFrom: ${return "N";}
      vep_cache: vep_cache
      output_basename: output_basename
      select_vars_mode: select_vars_mode
    out:
      [mutect2_filtered_stats, mutect2_filtered_vcf, mutect2_vep_vcf, mutect2_vep_tbi, mutect2_vep_maf]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4