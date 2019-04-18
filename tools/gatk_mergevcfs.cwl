cwlVersion: v1.0
class: CommandLineTool
id: gatk4_mergevcfs
label: GATK Merge VCF
doc: "Merge input vcfs"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 2
baseCommand: [/gatk, MergeVcfs]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx2000m"
      --TMP_DIR=./TMP
      --CREATE_INDEX=true
      --SEQUENCE_DICTIONARY=$(inputs.reference_dict.path)
      ${
        var cmd = "--OUTPUT=" + inputs.output_basename + "." + inputs.tool_name + ".merged.vcf.gz "
        if (typeof inputs.silent_flag !== 'undefined' && inputs.silent_flag == 1){
          cmd += "--VALIDATION_STRINGENCY SILENT"
        }
        return cmd
      }
      

inputs:
  input_vcfs:
    type:
      type: array
      items: File
      inputBinding:
        prefix: -I
    secondaryFiles: [.tbi]
    inputBinding:
      position: 1
  reference_dict: File
  tool_name: string
  output_basename: string
  silent_flag: ['null', int]
outputs:
  merged_vcf:
    type: File
    outputBinding:
      glob: '*.merged.vcf.gz'
    secondaryFiles: [.tbi]
