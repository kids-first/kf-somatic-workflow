cwlVersion: v1.2
class: Workflow
id: prepare_regions
doc: |
  Workflow for variable interval manipulation. This workflow will do the
  following:
  1. If calling_padding provided, pad the calling regions
  2. If supplement_vcfs provided, concat those records to 1
  3. If blacklist_regions provided, subtract those regions from 2
  4. If scatter_count > 0, scatter 3's regions into <scatter_count> files

  Any or all of these steps can be skipped. If the user only wishes to apply a
  blacklist and scatter, step 3 will be applied to the calling regions.

  The final return of the workflow are:
  - processed regions in BED and INTERVALLIST formats
  - if scattering, the scattered regions in BED and INTERVALLIST formats
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
  calling_padding: { type: 'int?', doc: "Padding to add to the calling intervals." }
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
  prescatter_bed: { type: 'File', outputSource: gatk_intervallisttobed_prescatter/outbed }
  prescatter_bedgz: { type: 'File', outputSource: gatk_intervallisttobed_prescatter/outbedgz }
  scattered_intervallists: { type: 'File[]?', outputSource: gatk_intervallisttools_scatter/output }
  scattered_beds: { type: 'File[]?', outputSource: gatk_intervallisttobed_scatter/outbed }

steps:
  gatk_bedtointervallist_calling:
    run: ../tools/gatk_bedtointervallist.cwl
    when: $(inputs.input_bed.basename.search(/(bed|bed.gz)$/) != -1)
    in:
      input_bed: calling_regions
      output_filename:
        valueFrom: $(inputs.input_bed.basename.replace(/(.bed|.bed.gz)$/,".interval_list"))
      reference_dict: reference_dict
    out: [output]

  gatk_bedtointervallist_blacklist:
    run: ../tools/gatk_bedtointervallist.cwl
    when: $(inputs.input_bed != null && inputs.input_bed.basename.search(/(bed|bed.gz)$/) != -1)
    in:
      input_bed: blacklist_regions
      output_filename:
        valueFrom: $(inputs.input_bed.basename.replace(/(.bed|.bed.gz)$/,".interval_list"))
      reference_dict: reference_dict
    out: [output]

  gatk_intervallisttools_pad_calling:
    run: ../tools/gatk_intervallisttools.cwl
    when: $(inputs.padding != null)
    in:
      input:
        source: [gatk_bedtointervallist_calling/output, calling_regions]
        valueFrom: |
          $([].concat(self.filter(function(e) { return e != null })[0]))
      padding: calling_padding
      output_name:
        source: calling_regions
        valueFrom: |
          $(self.basename.replace(/(.bed|.bed.gz|.interval_list)$/,".padded.interval_list"))
    out: [output, count_output]

  gatk_intervallisttools_concat_supplement:
    run: ../tools/gatk_intervallisttools.cwl
    when: $(inputs.input[3] != null)
    in:
      input:
        source: [gatk_intervallisttools_pad_calling/output, gatk_bedtointervallist_calling/output, calling_regions, supplement_vcfs]
        valueFrom: |
          $(self[3].concat(self.filter(function(e) { return e != null})[0]))
      unique:
        valueFrom: $(1 == 1)
      include_filtered:
        valueFrom: $(1 == 1)
      output_name:
        source: calling_regions
        valueFrom: $(self.basename.replace(/(.bed|.bed.gz|.interval_list)$/,".supplemented.interval_list"))
    out: [output, count_output]

  gatk_intervallisttools_subtract_blacklist:
    run: ../tools/gatk_intervallisttools.cwl
    when: $(inputs.second_input.some(function(e) { return e != null }))
    in:
      input:
        source: [gatk_intervallisttools_concat_supplement/output, gatk_intervallisttools_pad_calling/output, gatk_bedtointervallist_calling/output, calling_regions]
        valueFrom: |
          $([].concat(self.filter(function(e) { return e != null })[0]))
      action:
        valueFrom: "SUBTRACT"
      second_input:
        source: [gatk_bedtointervallist_blacklist/output, blacklist_regions]
        valueFrom: |
          $([].concat(self.filter(function(e) { return e != null })[0]))
      output_name:
        source: calling_regions
        valueFrom: $(self.basename.replace(/(.bed|.bed.gz|.interval_list)$/,".blacklisted.interval_list"))
    out: [output, count_output]

  clt_pass_interval:
    run: ../tools/clt_pass_file.cwl
    in:
      infile:
        source: [gatk_intervallisttools_subtract_blacklist/output, gatk_intervallisttools_concat_supplement/output, gatk_intervallisttools_pad_calling/output, gatk_bedtointervallist_calling/output, calling_regions]
        valueFrom: |
          $([].concat(self.filter(function(e) { return e != null })[0])[0])
    out: [outfile]

  gatk_intervallisttobed_prescatter:
    run: ../tools/gatk_intervallisttobed.cwl
    in:
      input_intervallist: clt_pass_interval/outfile
      output_filename:
        valueFrom: $(inputs.input_intervallist.nameroot).bed
    out: [outbed, outbedgz]

  gatk_intervallisttools_scatter:
    run: ../tools/gatk_intervallisttools.cwl
    when: $(inputs.scatter_count > 0)
    in:
      input:
        source: [gatk_intervallisttools_subtract_blacklist/output, gatk_intervallisttools_concat_supplement/output, gatk_intervallisttools_pad_calling/output, gatk_bedtointervallist_calling/output, calling_regions]
        valueFrom: |
          $([].concat(self.filter(function(e) { return e != null })[0]))
      break_bands_at_multiples_of: break_bands_at_multiples_of
      scatter_count: scatter_count
      subdivision_mode: subdivision_mode
      unique:
        valueFrom: $(1 == 1)
    out: [output, count_output]

  gatk_intervallisttobed_scatter:
    run: ../tools/gatk_intervallisttobed.cwl
    when: $(inputs.input_intervallist != null)
    scatter: [input_intervallist]
    in:
      input_intervallist: gatk_intervallisttools_scatter/output
      output_filename:
        valueFrom: $(inputs.input_intervallist.nameroot).bed
    out: [outbed, outbedgz]

$namespaces:
  sbg: https://sevenbridges.com
"sbg:license": Apache License 2.0
"sbg:publisher": KFDRC
hints:
  - class: 'sbg:AWSInstanceType'
    value: c5.2xlarge
