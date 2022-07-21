cwlVersion: v1.2
class: Workflow
id: kfdrc-amplicon-architect-workflow
label: Kids First DRC Amplicon Architect Workflow
doc: "ecDNA"

requirements:
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement
- class: ScatterFeatureRequirement

inputs:
  # Required
  aa_data_repo: { type: File, doc: "Reference tar ball obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/"}
  aa_data_ref_version: { type: ['null', {type: enum, name: wgs_mode, symbols: ["GRCh38", "hg19", "GRCh37", "mm10", "GRCm38"]}], doc: "Genome version in data repo to use", default: "GRCh38" }
  tumor_align_file: { type: File, doc: "tumor read alignment file. Can be cram if .cns file available", secondaryFiles: [ '^.bai?', '.crai?' ] }
  output_basename: { type: string, doc: "File name prefix for steps that require it" }
  mosek_license_file: { type: File, doc: "This tool uses some software that requires a license file. You can get a personal or institutional one from https://www.mosek.com/license/request/." }
  # If cram input and cnvkit not run
  ref_cache: { type: 'File?', doc: "For cram input, provide tar ball of output of running misc/seq_cache_populate.pl from samtools on the reference fasta" }
  reference: { type: 'File?', doc: "fasta file, needed if cram input or cnv kit cnn not already built", secondaryFiles: [.fai]}
  # If proper .cns file exists
  cnvkit_cns: { type: 'File?', doc: ".cns file from previous CNVkit run, if available. DO NOT USE .call.cns" }
  # If CNVkit needs running and .cnn file exists
  cnvkit_cnn: { type: 'File?', doc: "If running using an existing .cnn, supply here" }
  # If CNVkit needs running and .cnn NOT available
  normal_align_file: { type: 'File?', doc: "tumor read alignment file. Can be cram if .cns file available", secondaryFiles: [ { pattern: '^.bai', required: false },  { pattern: '.crai', required: false } ] }
  annotation_file: { type: 'File?', doc: "refFlat.txt file, needed if cnv kit cnn not already built" }
  male_input_flag: { type: 'boolean?', doc: "Is the input male?", default: false }

outputs:
  aa_cnv_seeds: { type: File,  doc: "Bed file with candidate regions to search", outputSource: prepare_aa/aa_cnv_seeds }
  aa_summary: { type: File, doc: "summary for all amplicons detected by AA", outputSource: amplicon_architect/summary }
  aa_cycles: { type: 'File[]', doc: "text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts", outputSource: amplicon_architect/cycles }
  graph: { type: 'File[]', doc: 'A text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts', outputSource: amplicon_architect/graph }
  sv_pdf: { type: 'File[]', doc: "PDF image file displaying the SV view of AA", outputSource: amplicon_architect/sv_pdf }
  sv_png: { type: 'File[]', doc: "PNG image file displaying the SV view of AA", outputSource: amplicon_architect/sv_png }
  amplicon_classification_profiles: { type: 'File[]?', doc: "abstract classification of the amplicon", outputSource: amplicon_classifier/amplicon_classification_profiles }
  gene_list: { type: 'File[]?', doc: "genes present on amplicons with each classification", outputSource: amplicon_classifier/gene_list }
  ecDNA_cts: { type: 'File[]?', doc: "not yet defined", outputSource:  amplicon_classifier/ecDNA_cts }

steps:
  samtools_tumor_cram_to_bam:
    run: ../tools/samtools_calmd.cwl
    when: $(inputs.input_reads.nameext == ".cram")
    in:
      input_reads: tumor_align_file
      reference: reference
    out: [bam_file]

  samtools_normal_cram_to_bam:
    run: ../tools/samtools_calmd.cwl
    when: $(inputs.input_reads != null && inputs.input_reads.nameext == ".cram")
    in:
      input_reads: normal_align_file
      reference: reference
    out: [bam_file]

  untar_data_repo:
    run: ../tools/untar_gzip.cwl
    in:
      input_tar: aa_data_repo
      out_dir_name:
        valueFrom: ${return "data_repo"}
    out: [output]


  cnvkit_batch:
    run: ../tools/cnvkit_batch_only.cwl
    when: $(inputs.cnvkit_cnn != null || inputs.input_control != null)
    in:
      input_sample:
        source: [samtools_tumor_cram_to_bam/bam_file, tumor_align_file]
        pickValue: first_non_null
      input_control:
        source: [samtools_normal_cram_to_bam/bam_file, normal_align_file]
        pickValue: first_non_null
      reference:
        # samtools needs fasta, can't use fasta if cnn file provided
        source: [cnvkit_cnn, reference]
        valueFrom: "${
          if (self[0] == null){
            return self[1];
          }
          else{
            return null;
          }
        }"
      cnvkit_cnn: cnvkit_cnn
      annotation_file: annotation_file
      male_input_flag: male_input_flag
      output_basename: output_basename
    out: [output_cnr, output_cns, output_scatter, output_diagram, output_cnn]

  cns_to_aa_bed:
    run: ../tools/cns_to_aa_bed.cwl
    in:
      input_cns:
        source: [cnvkit_batch/output_cns, cnvkit_cns]
        pickValue: first_non_null
    out: [aa_cns_bed]

  prepare_aa:
    run: ../tools/prepare_aa.cwl
    in:
      data_repo: untar_data_repo/output
      data_ref_version: aa_data_ref_version
      sample: output_basename
      sorted_bam:
        source: [samtools_tumor_cram_to_bam/bam_file, tumor_align_file]
        pickValue: first_non_null
      cnv_bed: cns_to_aa_bed/aa_cns_bed
      ref_cache: ref_cache
    out: [aa_cnv_seeds, coverage_stats]

  amplicon_architect:
    run: ../tools/amplicon_architect.cwl
    in:
      data_repo: untar_data_repo/output
      data_ref_version: aa_data_ref_version
      bam:
        source: [samtools_tumor_cram_to_bam/bam_file, tumor_align_file]
        pickValue: first_non_null
      output_basename: output_basename
      coverage_stats: prepare_aa/coverage_stats
      aa_bed: prepare_aa/aa_cnv_seeds
      ref_cache: ref_cache
      mosek_license_file: mosek_license_file
    out: [summary, cycles, graph, sv_pdf, sv_png]

  amplicon_classifier:
    run: ../tools/aa_classifier.cwl
    in:
      data_repo: untar_data_repo/output
      data_ref_version: aa_data_ref_version
      cycles: amplicon_architect/cycles
      graph: amplicon_architect/graph
    scatter: [cycles, graph]
    scatterMethod: dotproduct
    out: [amplicon_classification_profiles, gene_list, ecDNA_cts]

$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 3