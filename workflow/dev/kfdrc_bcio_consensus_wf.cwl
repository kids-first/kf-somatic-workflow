cwlVersion: v1.0
class: Workflow
id: kfdrc_bcbio_consensus_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  input_vcf:
    type:
        type: array
        items:
            type: array
            items: File
    doc: "Array of array of vcf files to consensus call"
  min_overlap:
    type: ['null', int]
    default: 2
    doc: "Min number of callers to declare consensus.  Default is 2"
  reference: File
  output_basename: string[]
  tool_name_csv: {type: string[], doc: "csv string with tools used.  should be in same order as input file array"}

outputs:
  consensus_vcf: {type: 'File[]', outputSource: bcbio_ensemble/consensus_vcf}

steps:
  bcbio_ensemble:
    run: ../../tools/bcbio_variant_recall_ensemble.cwl
    in:
      input_vcf: input_vcf
      reference: reference
      output_basename: output_basename
      tool_name_csv: tool_name_csv
      min_overlap: min_overlap
    scatter: 
        - input_vcf
        - output_basename
        - tool_name_csv
    scatterMethod: dotproduct
    out: [consensus_vcf]


$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:AWSInstanceType'
    value: c5.9xlarge;ebs-gp2;400
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4