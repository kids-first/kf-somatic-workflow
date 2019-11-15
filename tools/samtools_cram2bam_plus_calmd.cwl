cwlVersion: v1.0
class: CommandLineTool
id: samtools_cram2bam_plus_calmd
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 12000
    coresMin: $(inputs.threads)
  
baseCommand: ["/bin/bash -c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      ${
        var bam_name = inputs.input_reads.nameroot + ".bam";
        var cmd = "samtools view -@ " + inputs.threads + " -h -T " + inputs.reference.path + " " + inputs.input_reads.path
        + " | samtools calmd -@ " + inputs.threads + " -b --reference " + inputs.reference.path + " - > " + bam_name + ";";
        if(inputs.input_reads.basename == bam_name){
          cmd = ">&2 echo input reads already have bam extension, indexing and passing through; cp " + inputs.input_reads.path
          + " " + bam_name + ";"
        }
        cmd += "samtools index -@ " + inputs.threads + " " + bam_name + " " + inputs.input_reads.nameroot + ".bai;"
        return cmd;
      }
inputs:
  input_reads: File
  threads:
    type: ['null', int]
    default: 16
  reference: {type: File, secondaryFiles: [.fai]}
outputs:
  bam_file:
    type: File
    outputBinding:
      glob: '*.bam'
    secondaryFiles: [^.bai]
