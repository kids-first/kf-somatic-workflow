cwlVersion: v1.2
class: Workflow
id: kfdrc-amplicon-architect-workflow

label: Kids First DRC Amplicon Architect Workflow
doc: |
  # Kids First DRC Amplicon Architect(AA) Workflow
  This is an extra-chromosomal DNA (ecDNA) pipeline designed to detect ecDNA candidates.
  Briefly, chromosomal fragments in tumor cells break off to form circular DNA.
  This DNA acts sort of like a plasmid/mtDNA in that it is not beholden to the usual rules of cell division and can have massive copy numbers in each cell.
  The tools for this pipeline were derived from https://github.com/AmpliconSuite/AmpliconSuite-pipeline.

  ## [AmpliconSuite-pipeline](../tools/ampliconsuite_pipeline.cwl)
  The primary work element of the workflow is a CWL wrapper for [AmpliconSuite-pipeline.py](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/v1.5.2/AmpliconSuite-pipeline.py).
  The script is itself an intelligent pipeline that enables all steps (alignment, CNV calling, seed interval detection) prior to running AmpliconArchitect, and invokes AmpliconClassifier.
  This tool can reasonably be run as a standalone element to achieve all the same results that one would get from the workflow.

  ## [Workflow](../workflow/kfdrc_production_amplicon_architect.cwl)
  As mentioned above the AmpliconSuite-pipeline can simply be run to get results, so why a workflow? Mostly it comes down to minor timecost considerations. They are as follows:
  - While modern CNVkit can handle CRAMs, it takes considerably longer to process them than BAMs. Even when accounting for the time to convert the CRAMs to BAMs before CNVkit, processing CRAMs takes longer and costs 20-30% more.
  - While AmpliconSuite-pipeline can run both CNVkit and AmpliconArchitect in the same run, the pipeline breaks these up into separate steps. The reason for this is twofold:
    - CNVkit is multithreaded, AmpliconArchitect is single threaded. To limit the number of idle CPUs, we switch to a smaller instance for AmpliconArchitect
    - Both of these steps can be quite long and are vulnerable to spot instance kills. Breaking these up into separate steps reduces that vulnerability and better leverages memoization
  ### Inputs (critical):
   - `aa_data_repo`: Reference tar ball obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/
   - `tumor_align_file`: "Tumor read alignment file. Can be BAM or CRAM and must include index (BAI/CRAI)"
   - `output_basename`: "File name prefix for steps that require it"
   - `mosek_license_file`: "This workflow uses some software that requires a license file. You can get a personal or institutional one from https://www.mosek.com/license/request/."
  ### Input (if CRAM input):
   - `reference`: Reference fasta file
  ### Input (if `.cns` file available):
   - `cnvkit_cns`: .cns file from previous CNVkit run, if available. DO NOT USE .call.cns
  ### Input (if no CNVkit inputs available):
   - `normal_align_file`: Normal read alignment file. Can be CRAM or BAM
  ### Outputs:
   - `aa_cnv_seeds`: Bed file with candidate regions to search
   - `aa_summary`: Summary for all amplicons detected by AA
   - `aa_cycles`: Text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
   - `aa_graph`: A text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
   - `aa_sv_png`: PNG image file displaying the SV view of AA
   - `aa_classification_profiles`: Abstract classification of the amplicon
   - `aa_gene_list`: Genes present on amplicons with each classification
requirements:
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
  # Required
  aa_data_repo: {type: File, doc: "Reference tar ball obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/"}
  tumor_align_file: {type: File, doc: "Tumor read alignment file. Can be cram or bam", secondaryFiles: [{pattern: '^.bai', 
          required: false}, {pattern: '.crai', required: false}]}
  samtools_calmd_threads: {type: 'int?', default: 8}
  output_basename: {type: string, doc: "File name prefix for steps that require it"}
  mosek_license_file: {type: File, doc: "This tool uses some software that requires a license file. You can get a personal or institutional
      one from https://www.mosek.com/license/request/."}
  # If cram input and cnvkit not run
  reference: {type: 'File?', doc: "fasta file, needed if cram input or cnv kit cnn not already built", secondaryFiles: [.fai]}
  # If proper .cns file exists
  cnvkit_cns: {type: 'File?', doc: ".cns file from previous CNVkit run, if available. DO NOT USE .call.cns"}
  # If CNVkit needs running and .cnn NOT available
  normal_align_file: {type: 'File?', doc: "Normal read alignment file. Can cram or bam", secondaryFiles: [{pattern: '^.bai', 
          required: false}, {pattern: '.crai', required: false}]}
  cnvkit_cpu: {type: 'int?', doc: "CPUs to allocate to CNVkit running"}
  cnvkit_ram: {type: 'int?', doc: "GB of RAM to allocate to CNVkit running"}

outputs:
  aa_cnv_seeds: {type: 'File?', doc: "Bed file with candidate regions to search", outputSource: ampliconsuite_cnvkit/seeds}
  aa_summary: {type: File, doc: "summary for all amplicons detected by AA", outputSource: ampliconsuite_aa_ac/aa_summary}
  aa_cycles: {type: 'File[]', doc: "text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence,
      discordant, concordant, source) and their copy counts", outputSource: ampliconsuite_aa_ac/aa_cycles}
  aa_graph: {type: 'File[]', doc: 'A text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence,
      discordant, concordant, source) and their copy counts', outputSource: ampliconsuite_aa_ac/aa_graph}
  aa_sv_png: {type: 'File[]', doc: "PNG image file displaying the SV view of AA", outputSource: ampliconsuite_aa_ac/aa_sv_png}
  aa_classification_profiles: {type: 'File[]?', doc: "abstract classification of the amplicon", outputSource: 
      ampliconsuite_aa_ac/ac_profiles}
  aa_gene_list: {type: 'File[]?', doc: "genes present on amplicons with each classification", outputSource: 
      ampliconsuite_aa_ac/ac_gene_list}

steps:
  samtools_calmd_tumor_cram_to_bam:
    hints:
    - class: 'sbg:AWSInstanceType'
      value: c6i.2xlarge
    run: ../tools/samtools_calmd.cwl
    when: $(inputs.input_reads.nameext == ".cram" && inputs.normal_reads != null)
    in:
      normal_reads:
        source: normal_align_file
        valueFrom: |
          $(self ? 1 : null)
      input_reads: tumor_align_file
      reference: reference
      threads: samtools_calmd_threads
    out: [bam_file]
  samtools_calmd_normal_cram_to_bam:
    hints:
    - class: 'sbg:AWSInstanceType'
      value: c6i.2xlarge
    run: ../tools/samtools_calmd.cwl
    when: $(inputs.input_reads != null && inputs.input_reads.nameext == ".cram")
    in:
      input_reads: normal_align_file
      reference: reference
      threads: samtools_calmd_threads
    out: [bam_file]
  ampliconsuite_cnvkit:
    run: ../tools/ampliconsuite_pipeline.cwl
    when: $(inputs.normal_bam != null && inputs.cnv_bed == null)
    in:
      data_repo: aa_data_repo
      mosek_license_file: mosek_license_file
      sorted_tumor_bam:
        source: [samtools_calmd_tumor_cram_to_bam/bam_file, tumor_align_file]
        pickValue: first_non_null
      normal_bam:
        source: [samtools_calmd_normal_cram_to_bam/bam_file, normal_align_file]
        pickValue: first_non_null
      cnv_bed: cnvkit_cns
      sample_name: output_basename
      extra_args:
        valueFrom: >-
          --no_QC
      cpu: cnvkit_cpu
      ram: cnvkit_ram
    out: [log, timing_log, run_meta, sample_meta, seeds, cnvkit_results, aa_results, ac_results, aa_summary, aa_graph, aa_cycles,
      aa_sv_png, aa_sv_pdf, ac_profiles, ac_gene_list]
  ampliconsuite_aa_ac:
    hints:
    - class: 'sbg:AWSInstanceType'
      value: m6i.large
    run: ../tools/ampliconsuite_pipeline.cwl
    in:
      data_repo: aa_data_repo
      mosek_license_file: mosek_license_file
      sorted_tumor_bam:
        source: [samtools_calmd_tumor_cram_to_bam/bam_file, tumor_align_file]
        pickValue: first_non_null
      cram_reference: reference
      cnv_bed:
        source: [ampliconsuite_cnvkit/seeds, cnvkit_cns]
        pickValue: first_non_null
      sample_name: output_basename
      extra_args:
        valueFrom: >-
          --run_AA --run_AC
      cpu:
        valueFrom: $(1)
      ram:
        valueFrom: $(8)
    out: [log, timing_log, run_meta, sample_meta, seeds, cnvkit_results, aa_results, ac_results, aa_summary, aa_graph, aa_cycles,
      aa_sv_png, aa_sv_pdf, ac_profiles, ac_gene_list]
$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 2
"sbg:license": Apache License 2.0
"sbg:publisher": KFDRC
