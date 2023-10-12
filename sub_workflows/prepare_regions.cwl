cwlVersion: v1.2
class: Workflow
id: prepare_regions
doc: |
  Given a set of calling regions, in BED or INTERVALLIST format, blacklist
  regions, in BED or INTERVALLIST format, and scatter parameters, this workflow
  will subtract the blacklist regions from the calling reigons and return a list
  of roughly equal sized BEDs and INTERVALLISTs containing the remaining regions.
requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
inputs:
  calling_regions: { type: 'File', doc: "BED or INTERVALLIST file containing a set of genomic regions." }
  supplement_vcfs: { type: 'File[]?', secondaryFiles: [{ pattern: '.tbi', required: false }, { pattern: '.csi', required: false }], doc: "Supplemental variant calls in BCF/VCF(.GZ) format with associated index to add to the calling regions." }
  blacklist_regions: { type: 'File?', doc: "BED or INTERVALLIST file containing a set of genomic regions to remove from the calling regions." }
  reference_dict: { type: 'File?', doc: "Reference sequence dictionary, or BAM/VCF/IntervalList from which a dictionary can be extracted. Required if calling_regions are in a BED file." }
  break_bands_at_multiples_of: { type: 'int?', default: 1000000, doc: "If set to a positive value will create a new interval list with the original intervals broken up at integer multiples of this value. Set to 0 to NOT break up intervals." }
  scatter_count: { type: 'int?', default: 50, doc: "Total number of scatter intervals and beds to make" }
  subdivision_mode:
    type:
      - 'null'
      - type: enum
        name: subdivision_mode
        symbols: [ "INTERVAL_SUBDIVISION", "BALANCING_WITHOUT_INTERVAL_SUBDIVISION", "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW", "INTERVAL_COUNT", "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER" ]
    default: "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
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
outputs:
  prescatter_intervallist: { type: 'File', outputSource: clt_pass_interval/outfile }
  prescatter_bed: { type: 'File', outputSource: gatk_intervallisttobed_prescatter/output }
  scattered_intervallists: { type: 'File[]', outputSource: gatk_intervallisttools_scatter/output }
  scattered_beds: { type: 'File[]', outputSource: gatk_intervallisttobed_scatter/output }

steps:
  gatk_bedtointervallist_calling:
    run: ../tools/gatk_bedtointervallist.cwl
    when: $(inputs.input_bed.basename.search(/(bed|bed.gz)$/) != -1)
    in:
      input_bed: calling_regions
      output_filename:
        valueFrom: $(inputs.input_bed.nameroot).interval_list
      reference_dict: reference_dict
    out: [output]

  gatk_bedtointervallist_blacklist:
    run: ../tools/gatk_bedtointervallist.cwl
    when: $(inputs.input_bed != null && inputs.input_bed.basename.search(/(bed|bed.gz)$/) != -1)
    in:
      input_bed: blacklist_regions
      output_filename:
        valueFrom: $(inputs.input_bed.nameroot).interval_list
      reference_dict: reference_dict
    out: [output]

  gatk_intervallisttools_concat_supplement:
    run: ../tools/gatk_intervallisttools.cwl
    when: $(inputs.input[2] != null)
    in:
      input:
        source: [gatk_bedtointervallist_calling/output, calling_regions, supplement_vcfs]
        valueFrom: |
          $(self[0] != null ? self[2].concat(self[0]) : self[2].concat(self[1]))
      unique:
        valueFrom: $(1 == 1)
      output_name:
        source: calling_regions
        valueFrom: $(self.basename.replace(/(.bed|.bed.gz|.interval_list)$/,".combined.interval_list"))
    out: [output, count_output]

  gatk_intervallisttools_subtract_blacklist:
    run: ../tools/gatk_intervallisttools.cwl
    when: $(inputs.second_input.some(function(e) { return e != null }))
    in:
      input:
        source: [gatk_intervallisttools_concat_supplement/output, gatk_bedtointervallist_calling/output, calling_regions]
        valueFrom: |
          $(self[0] != null ? self[0] : self[1] != null ? [self[1]] : [self[2]])
      action:
        valueFrom: "SUBTRACT"
      second_input:
        source: [gatk_bedtointervallist_blacklist/output, blacklist_regions]
        valueFrom: |
          $(self[0] != null ? [self[0]] : [self[1]])
      output_name:
        source: calling_regions
        valueFrom: $(self.basename.replace(/(.bed|.bed.gz|.interval_list)$/,".blacklisted.interval_list"))
    out: [output, count_output]

  clt_pass_interval:
    run: ../tools/clt_pass_file.cwl
    in:
      infile:
        source: [gatk_intervallisttools_subtract_blacklist/output, gatk_intervallisttools_concat_supplement/output, gatk_bedtointervallist_calling/output, calling_regions]
        valueFrom: |
          $(self[0] != null ? self[0][0] : self[1] != null ? self[1][0] : self[2] != null ? self[2]: self[3])
    out: [outfile]

  gatk_intervallisttobed_prescatter:
    run: ../tools/gatk_intervallisttobed.cwl
    in:
      input_intervallist: clt_pass_interval/outfile
      output_filename:
        valueFrom: $(inputs.input_intervallist.nameroot).bed
    out: [output]

  gatk_intervallisttools_scatter:
    run: ../tools/gatk_intervallisttools.cwl
    in:
      input:
        source: [gatk_intervallisttools_subtract_blacklist/output, gatk_intervallisttools_concat_supplement/output, gatk_bedtointervallist_calling/output, calling_regions]
        valueFrom: |
          $(self[0] != null ? self[0] : self[1] != null ? self[1] : self[2] != null ? [self[2]] : [self[3]])
      break_bands_at_multiples_of: break_bands_at_multiples_of
      scatter_count: scatter_count
      subdivision_mode: subdivision_mode
      unique:
        valueFrom: $(1 == 1)
    out: [output, count_output]

  gatk_intervallisttobed_scatter:
    run: ../tools/gatk_intervallisttobed.cwl
    scatter: [input_intervallist]
    in:
      input_intervallist: gatk_intervallisttools_scatter/output
      output_filename:
        valueFrom: $(inputs.input_intervallist.nameroot).bed
    out: [output]

$namespaces:
  sbg: https://sevenbridges.com
"sbg:license": Apache License 2.0
"sbg:publisher": KFDRC
hints:
  - class: 'sbg:AWSInstanceType'
    value: c5.2xlarge
