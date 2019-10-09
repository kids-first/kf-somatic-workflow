cwlVersion: v1.0
class: CommandLineTool
id: strelka2
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 10000
    coresMin: 36
  - class: DockerRequirement
    dockerPull: 'obenauflab/strelka'

baseCommand: [/strelka-2.9.3.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --normalBam $(inputs.input_normal_aligned.path)
      --tumorBam $(inputs.input_tumor_aligned.path)
      --ref $(inputs.reference.path)
      --callRegions $(inputs.hg38_strelka_bed.path)
      ${
        var arg = "--runDir=./";
        if (inputs.exome_flag == 'Y'){
          arg += " --exome"
        }
        return arg
      } && ./runWorkflow.py
      -m local
      -j 36

inputs:
  reference: {type: File, secondaryFiles: [.fai]}
  hg38_strelka_bed: { type: File, secondaryFiles: [.tbi], label: gzipped bed file }
  exome_flag: { type: ['null', string], doc: "Y if exome/capture, defaults to WGS"}
  input_tumor_aligned:
    type: File
    secondaryFiles: |
      ${
        var dpath = self.location.replace(self.basename, "")
        if(self.nameext == '.bam'){
          return {"location": dpath+self.nameroot+".bai", "class": "File"}
        }
        else{
          return {"location": dpath+self.basename+".crai", "class": "File"}
        }
      }
    doc: "tumor BAM or CRAM"
  input_normal_aligned:
    type: File
    secondaryFiles: |
      ${
        var dpath = self.location.replace(self.basename, "")
        if(self.nameext == '.bam'){
          return {"location": dpath+self.nameroot+".bai", "class": "File"}
        }
        else{
          return {"location": dpath+self.basename+".crai", "class": "File"}
        }
      }
    doc: "normal BAM or CRAM"
outputs:
  output_snv:
    type: File
    outputBinding:
      glob: 'results/variants/*.snvs.vcf.gz'
    secondaryFiles: [.tbi]
  output_indel:
    type: File
    outputBinding:
      glob: 'results/variants/*.indels.vcf.gz'
    secondaryFiles: [.tbi]

