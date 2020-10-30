cwlVersion: v1.0
class: CommandLineTool
id: cnvkit_export_theta2
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: 4

baseCommand: [cnvkit.py, export, theta]  
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >- 
      -r $(inputs.reference_cnn.path)
      -v $(inputs.paired_vcf.path)
      -i $(inputs.tumor_ID)
      -n $(inputs.normal_ID)
      $(inputs.tumor_cns.path)

inputs:
  tumor_cns: File
  reference_cnn: File
  paired_vcf: File
  normal_ID: string
  tumor_ID: string

outputs:
  call_interval_count:
    type: File
    outputBinding:
      glob: '*.call.interval_count'
  call_tumor_snp:
    type: File
    outputBinding:
      glob: '*.call.tumor.snp_formatted.txt'
  call_normal_snp:
    type: File
    outputBinding:
      glob: '*.call.normal.snp_formatted.txt'
    
