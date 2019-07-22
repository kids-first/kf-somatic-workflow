cwlVersion: v1.0
class: CommandLineTool
id: cnvkit
requirements: 
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement 
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3'
  - class: ResourceRequirement
    ramMin: 32000
    coresMin: 16
arguments: 
  - position: 1
    shellQuote: false
    valueFrom: >-
      cnvkit.py
