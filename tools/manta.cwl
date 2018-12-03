cwlVersion: v1.0
class: CommandLineTool
id: kf-manta-sv
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 10000
    coresMin: 8
  - class: DockerRequirement
    dockerPull: 'kfdrc/manta:latest'

baseCommand: [/manta-1.4.0.centos6_x86_64/bin/configManta.py]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      ${
        if (typeof inputs.input_tumor_cram === 'undefined' || inputs.input_tumor_cram === null){
          return "--bam ".concat(inputs.input_normal_cram.path);
        }
        else if (typeof inputs.input_normal_cram === 'undefined' || inputs.input_normal_cram === null){
          return "--tumorBam ".concat(inputs.input_tumor_cram.path)
        }
        else{
          return "--tumorBam ".concat(inputs.input_tumor_cram.path, " --normalBam ", inputs.input_normal_cram.path)
        }
      }
      --ref $(inputs.reference.path)
      --callRegions $(inputs.ref_bed.path)
      --runDir=./ && ./runWorkflow.py
      -m local
      -j 8
      --quiet
      && mv results/variants/somaticSV.vcf.gz $(inputs.output_basename).somaticSV.vcf.gz
      && mv results/variants/somaticSV.vcf.gz.tbi $(inputs.output_basename).somaticSV.vcf.gz.tbi

inputs:
    reference: {type: File, secondaryFiles: [^.dict, .fai]}
    ref_bed: {type: File, secondaryFiles: [.tbi]}
    input_tumor_cram: {type: ["null", File], secondaryFiles: [.crai]}
    input_normal_cram: {type: ["null", File], secondaryFiles: [.crai]}
    output_basename: string
outputs:
  - id: output_sv
    type: File
    outputBinding:
      glob: '*.somaticSV.vcf.gz'
    secondaryFiles: [.tbi]
