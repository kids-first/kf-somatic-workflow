cwlVersion: v1.2
class: CommandLineTool
id: gatk_intervallisttools
doc: |
  This tool offers multiple interval list file manipulation capabilities,
  including: sorting, merging, subtracting, padding, and other set-theoretic
  operations. The default action is to merge and sort the intervals provided in
  the INPUTs. Other options, e.g. interval subtraction, are controlled by the
  arguments.
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.4.0.0'
baseCommand: [gatk, IntervalListTools]
arguments:
  - position: 10
    prefix: '&&'
    shellQuote: false
    valueFrom: >
      $(inputs.subdivision_mode != null ? "" : "exit 0;")
      find ./temp_* -name *interval_list | sed 'p;s#temp_\([0-9]\+\).*#\1_scattered.interval_list#' | xargs -P $(inputs.cpu) -n 2 mv
inputs:
  # INPUT ARGUMENTS
  input:
    type:
      type: array
      items: File
      inputBinding:
        prefix: --INPUT
    inputBinding:
      position: 2
    doc: |
      One or more interval lists. If multiple interval lists are provided the output
      is theresult of merging the inputs.
  action:
    type:
      - 'null'
      - type: enum
        name: action
        symbols: [ "CONCAT", "UNION", "INTERSECT", "SUBTRACT", "SYMDIFF", "OVERLAPS"]
    inputBinding:
      prefix: --ACTION
      position: 2
    doc: |
      Action to take on inputs.  Default value: CONCAT. CONCAT (The concatenation of
      all the intervals in all the INPUTs, no sorting or merging of
      overlapping/abutting intervals implied. Will result in a possibly unsorted list
      unless requested otherwise.) UNION (Like CONCATENATE but with UNIQUE and SORT
      implied, the result being the set-wise union of all INPUTS, with overlapping
      and abutting intervals merged into one.) INTERSECT (The sorted and merged set
      of all loci that are contained in all of the INPUTs.) SUBTRACT (Subtracts the
      intervals in SECOND_INPUT from those in INPUT. The resulting loci are those in
      INPUT that are not in SECOND_INPUT.) SYMDIFF (Results in loci that are in INPUT
      or SECOND_INPUT but are not in both.) OVERLAPS (Outputs the entire intervals
      from INPUT that have bases which overlap any interval from SECOND_INPUT. Note
      that this is different than INTERSECT in that each original interval is either
      emitted in its entirety, or not at all.)
  second_input:
    type:
      - 'null'
      -  type: array
         items: File
         inputBinding:
           prefix: --SECOND_INPUT
    inputBinding:
      position: 2
    doc: |
      Second set of intervals for SUBTRACT and DIFFERENCE operations. This argument
      may be specified 0 or more times.
  break_bands_at_multiples_of: { type: 'int?', inputBinding: { position: 2, prefix: "--BREAK_BANDS_AT_MULTIPLES_OF" }, doc: "If set to a positive value will create a new interval list with the original intervals broken up at integer multiples of this value. Set to 0 to NOT break up intervals." }
  padding: { type: 'int?', inputBinding: { position: 2, prefix: "--PADDING" }, doc: "The amount to pad each end of the intervals by before other operations are undertaken.  Negative numbers are allowed and indicate intervals should be shrunk. Resulting intervals < 0 bases long will be removed. Padding is applied to the interval lists (both INPUT and SECOND_INPUT, if provided) before the ACTION is performed." }

  # OUTPUT ARGUMENTS
  output_name: { type: 'string?', default: "out.interval_list", inputBinding: { position: 2, prefix: "--OUTPUT", valueFrom: '$(inputs.subdivision_mode != null ? "." : self)' }, doc: "The output interval list file to write (if SCATTER_COUNT == 1) or the directory into which to write the scattered interval sub-directories (if SCATTER_COUNT > 1)." }
  count_output_name: { type: 'string?', inputBinding: { position: 2, prefix: "--COUNT_OUTPUT" }, doc: "File to which to print count of bases or intervals in final output interval list.  When not set, value indicated by OUTPUT_VALUE will be printed to stdout.  If this parameter is set, OUTPUT_VALUE must not be NONE." }
  output_value:
    type:
      - 'null'
      - type: enum
        name: output_value
        symbols: [ "NONE", "BASES", "INTERVALS" ]
    inputBinding:
      prefix: --OUTPUT_VALUE
      position: 2
    doc: |
      What value to output to COUNT_OUTPUT file or stdout (for scripting).  If
      COUNT_OUTPUT is provided, this parameter must not be NONE.

  # SCATTERING ARGUMENTS
  scatter_content: { type: 'int?', inputBinding: { position: 2, prefix: "--SCATTER_CONTENT" }, doc: "When scattering with this argument, each of the resultant files will (ideally) have this amount of 'content', which  means either base-counts or interval-counts depending on SUBDIVISION_MODE. When provided, overrides SCATTER_COUNT." }
  scatter_count: { type: 'int?', default: 1, inputBinding: { position: 2, prefix: "--SCATTER_COUNT" }, doc: "The number of files into which to scatter the resulting list by locus; in some situations, fewer intervals may be emitted." }
  subdivision_mode:
    type:
      - 'null'
      - type: enum
        name: subdivision_mode
        symbols: [ "INTERVAL_SUBDIVISION", "BALANCING_WITHOUT_INTERVAL_SUBDIVISION", "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW", "INTERVAL_COUNT", "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER" ]
    inputBinding:
      prefix: --SUBDIVISION_MODE
      position: 2
    doc: |
      The mode used to scatter the interval list:
      - INTERVAL_SUBDIVISION (Scatter the interval list into similarly sized interval
        lists (by base count), breaking up intervals as needed.)
      - BALANCING_WITHOUT_INTERVAL_SUBDIVISION (Scatter the interval list into
        similarly sized interval lists (by base count), but without breaking up
        intervals.)
      - BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW (Scatter the interval
        list into similarly sized interval lists (by base count), but without
        breaking up intervals. Will overflow current interval list so that the
        remaining lists will not have too many bases to deal with.)
      - INTERVAL_COUNT (Scatter the interval list into similarly sized interval lists
        (by interval count, not by base count). Resulting interval lists will contain
        the same number of intervals except for the last, which contains the
        remainder.)
      - INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER (Scatter the interval list into
        similarly sized interval lists (by interval count, not by base count).
        Resulting interval lists will contain similar number of intervals.)

  include_filtered: { type: 'boolean?', default: false, inputBinding: { position: 2, prefix: "--INCLUDE_FILTERED", valueFrom: '$(self ? "true" : "false")' }, doc: "Whether to include filtered variants in the vcf when generating an interval list from vcf." }
  invert: { type: 'boolean?', inputBinding: { position: 2, prefix: "--INVERT", valueFrom: '$(self ? "true" : "false")' }, doc: "Produce the inverse list of intervals, that is, the regions in the genome that are <br>not</br> covered by any of the input intervals. Will merge abutting intervals first.  Output will be sorted." }
  sort: { type: 'boolean?', default: true, inputBinding: { position: 2, prefix: "--SORT", valueFrom: '$(self ? "true" : "false")' }, doc: "If true, sort the output interval list before writing it." }
  unique: { type: 'boolean?', default: false, inputBinding: { position: 2, prefix: "--UNIQUE", valueFrom: '$(self ? "true" : "false")' }, doc: "If true, unique the output interval list by merging overlapping regions, before writing it (implies sort=true)." }
  max_memory: { type: 'int?', default: 4, doc: "GB of RAM to allocate to the task." }
  cpu: { type: 'int?', default: 2, doc: "Minimum reserved number of CPU cores for the task." }
outputs:
  output: { type: 'File[]', outputBinding: { glob: '$(inputs.subdivision_mode != null ? "*.interval_list" : inputs.output_name)' } }
  count_output: { type: 'File?', outputBinding: { glob: $(inputs.count_output_name) } }
