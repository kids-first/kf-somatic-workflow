cwlVersion: v1.0
class: CommandLineTool
id: gatk_funcotatesegments
doc: "Functional annotation for segment files."
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory*1000)
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.4.1'
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      mkdir datasources_dir &&
      DATA_SOURCES_FOLDER="$PWD/datasources_dir" &&
      $( inputs.run_funcotatesegments ? 'tar xvzf '+inputs.data_sources_tgz.path+' -C datasources_dir --strip-components 1 && gatk' : 'echo gatk' )
  - position: 1
    shellQuote: false
    prefix: "--java-options"
    valueFrom: >-
      $("\"-Xmx"+Math.floor(inputs.max_memory*1000/1.074 - 1)+"M\"")
  - position: 2
    shellQuote: false
    valueFrom: >-
      FuncotateSegments
  - position: 3
    shellQuote: false
    prefix: "--data-sources-path"
    valueFrom: >-
      $DATA_SOURCES_FOLDER
  - position: 3
    shellQuote: false
    prefix: "-O"
    valueFrom: >-
      ${var pre = inputs.output_prefix ? inputs.output_prefix : inputs.segments ? inputs.segments.nameroot : 'output'; var ext = 'funcotated.tsv'; return pre+'.'+ext}
inputs:
  run_funcotatesegments:
    type: 'boolean?'
    doc: "If not true, this tool will produce no outputs."
  add_output_vcf_command_line:
    type: 'boolean?'
    doc: "If true, adds a command line header line to created VCF files."
    inputBinding:
      position: 3
      prefix: "--add-output-vcf-command-line"
  annotation_default:
    type:
      - 'null'
      - type: array
        items: string
        inputBinding:
          prefix: "--annotation-default"
    doc: "Annotations to include in all annotated variants if the annotation is not specified in the data sources (in the format :). This will add the specified annotation to every annotated variant if it is not already present."
    inputBinding:
      position: 3
  annotation_override:
    type:
      - 'null'
      - type: array
        items: string
        inputBinding:
          prefix: "--annotation-override"
    doc: "Override values for annotations (in the format :). Replaces existing annotations of the given name with given values."
    inputBinding:
      position: 3
  create_output_variant_index:
    type: 'boolean?'
    doc: "If true, create a VCF index when writing a coordinate-sorted VCF file."
    inputBinding:
      position: 3
      prefix: "--create-output-variant-index"
  create_output_variant_md5:
    type: 'boolean?'
    doc: "If true, create a a MD5 digest any VCF file created."
    inputBinding:
      position: 3
      prefix: "--create-output-variant-md5"
  data_sources_tgz:
    type: 'File?'
    doc: "TGZ File containing prebuild data sources folder."
  disable_bam_index_caching:
    type: 'boolean?'
    doc: "If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified. Caching is automatically disabled if there are no intervals specified."
    inputBinding:
      position: 3
      prefix: "--disable-bam-index-caching"
  disable_read_filter:
    type: 'boolean?'
    doc: "Read filters to be disabled before analysis."
    inputBinding:
      position: 3
      prefix: "--disable-read-filter"
  disable_sequence_dictionary_validation:
    type: 'boolean?'
    doc: "If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!"
    inputBinding:
      position: 3
      prefix: "--disable-sequence-dictionary-validation"
  disable_tool_default_read_filters:
    type: 'boolean?'
    doc: "Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)"
    inputBinding:
      position: 3
      prefix: "--disable-tool-default-read-filters"
  exclude_intervals:
    type:
      - 'null'
      - type: array
        items: string
        inputBinding:
          prefix: '-XL'
    inputBinding:
      position: 3
  exclude_interval_list:
    type: 'File?'
    doc: "One or more genomic intervals to exclude from processing. Use this input when providing interval list files or other file based inputs."
    inputBinding:
      position: 3
      prefix: "-XL"
  exclude_field:
    type:
      - 'null'
      - type: array
        items: string
        inputBinding:
          prefix: "--exclude-field"
    doc: "Fields that should not be rendered in the final output. Only exact name matches will be excluded."
    inputBinding:
      position: 3
  interval_exclusion_padding:
    type: 'int?'
    doc: "Amount of padding (in bp) to add to each interval you are excluding."
    inputBinding:
      position: 3
      prefix: "-ixp"
  interval_merging_rule:
    type:
      - 'null'
      - type: enum
        name: interval_merging_rule
        symbols: ["ALL","OVERLAPPING_ONLY"]
    doc: "By default, the program merges abutting intervals (i.e. intervals that are directly side-by-side but do not actually overlap) into a single continuous interval. However you can change this behavior if you want them to be treated as separate intervals instead."
    inputBinding:
      position: 3
      prefix: "-imr"
  interval_padding:
    type: 'int?'
    doc: "Amount of padding (in bp) to add to each interval you are including."
    inputBinding:
      position: 3
      prefix: "-ip"
  interval_set_rule:
    type:
      - 'null'
      - type: enum
        name: interval_set_rule
        symbols: ["UNION","INTERSECTION"]
    doc: "By default, the program will take the UNION of all intervals specified using -L and/or -XL. However, you can change this setting for -L, for example if you want to take the INTERSECTION of the sets instead. E.g. to perform the analysis only on chromosome 1 exomes, you could specify -L exomes.intervals -L 1 --interval-set-rule INTERSECTION. However, it is not possible to modify the merging approach for intervals passed using -XL (they will always be merged using UNION). Note that if you specify both -L and -XL, the -XL interval set will be subtracted from the -L interval set."
    inputBinding:
      position: 3
      prefix: "-isr"
  intervals:
    type:
      - 'null'
      - type: array
        items: string
        inputBinding:
          prefix: '-L'
    doc: "One or more genomic intervals over which to operate. Use this input when providing the intervals as strings."
    inputBinding:
      position: 3
  interval_list:
    type: 'File?'
    doc: "One or more genomic intervals over which to operate. Use this input when providing interval list files or other file based inputs."
    inputBinding:
      position: 3
      prefix: "-L"
  lenient:
    type: 'boolean?'
    doc: "Lenient processing of VCF files"
    inputBinding:
      position: 3
      prefix: "--lenient"
  lookahead_cache_bp:
    type: 'int?'
    doc: "Number of base-pairs to cache when querying variants. Can be overridden in individual data source configuration files."
    inputBinding:
      position: 3
      prefix: "--lookahead-cache-bp"
  output_file_format:
    type:
      - 'null'
      - type: enum
        name: output_file_format
        symbols: ["MAF","SEG","VCF"]
    default: "SEG"
    doc: "The output file format. Either VCF or MAF. Please note that MAF output for germline use case VCFs is unsupported."
    inputBinding:
      position: 3
      prefix: "--output-file-format"
  output_prefix:
    type: 'string?'
    doc: "String to use as the prefix for the outputs."
  read_filter:
    type: 'string[]?'
    doc: "Read filters to be applied before analysis"
    inputBinding:
      position: 3
      prefix: "--read-filter"
  read_index:
    type:
      - 'null'
      - type: array
        items: File
        inputBinding:
          prefix: "--read-index"
    doc: "Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically."
    inputBinding:
      position: 3
  read_validation_stringency:
    type:
      - 'null'
      - type: enum
        name: read_validation_stringency
        symbols: ["SILENT","STRICT","LENIENT"]
    doc: "Validation stringency for all SAM/BAM/CRAM/SRA files read by this program. The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded."
    inputBinding:
      position: 3
      prefix: "--read-validation-stringency"
  ref_version:
    type:
      - 'null'
      - type: enum
        name: ref_version
        symbols: ["hg19","hg38","b37"]
    doc: "The version of the Human Genome reference to use (e.g. hg19, hg38, etc.). This will correspond to a sub-folder of each data source corresponding to that data source for the given reference."
    inputBinding:
      position: 3
      prefix: "--ref-version"
  reference:
    type: 'File?'
    doc: "Reference fasta"
    secondaryFiles: ['.fai']
    inputBinding:
      position: 3
      prefix: "-R"
  remove_filtered_variants:
    type: 'boolean?'
    doc: "Ignore/drop variants that have been filtered in the input. These variants will not appear in the output file."
    inputBinding:
      position: 3
      prefix: "--remove-filtered-variants"
  segments:
    type: 'File?'
    doc: "Input segment file (tab-separated values). Must have a call column."
    inputBinding:
      position: 3
      prefix: "--segments"
  sequence_dictionary:
    type: 'File?'
    doc: "Use the given sequence dictionary as the master/canonical sequence dictionary. Must be a .dict file."
    inputBinding:
      position: 3
      prefix: "--sequence-dictionary"
  sites_only_vcf_output:
    type: 'boolean?'
    doc: "If true, don't emit genotype fields when writing vcf file output."
    inputBinding:
      position: 3
      prefix: "--sites-only-vcf-output"
  transcript_list:
    type:
      - 'null'
      - type: array
        items: string
        inputBinding:
          prefix: "--transcript-list"
    doc: "A set of transcript IDs to use for annotation to override selected transcript."
    inputBinding:
      position: 3
  transcript_list_file:
    type: 'File?'
    doc: "File to use as a list of transcripts (one transcript ID per line, version numbers are ignored)"
    inputBinding:
      position: 3
      prefix: "--transcript-list"
  transcript_selection_mode:
    type:
      - 'null'
      - type: enum
        name: transcript_selection_mode
        symbols: ["ALL","BEST_EFFECT","CANONICAL"]
    doc: "Method of detailed transcript selection. This will select the transcript for detailed annotation (CANONICAL, ALL, or BEST_EFFECT)."
    inputBinding:
      position: 3
      prefix: "--transcript-selection-mode"
  variant:
    type: 'File?'
    doc: "A VCF file containing variants"
    inputBinding:
      position: 3
      prefix: "--variant"
  minimum_segment_size:
    type: 'int?'
    default: 150
    doc: "Advanced Argument: The minimum number of bases for a variant to be annotated as a segment. Recommended to be changed only for use with FuncotateSegments."
    inputBinding:
      position: 3
      prefix: "--min-num-bases-for-segment-funcotation"
  max_memory:
    type: 'int?'
    default: 20
    doc: "Maximum GB of RAM to allocate for this tool."
outputs:
  funcotated_seg_simple_tsv: { type: 'File?', outputBinding: { glob: "*.funcotated.tsv" } }
  funcotated_gene_list_tsv: { type: 'File?', outputBinding: { glob: "*.funcotated.tsv.gene_list.txt" } }
