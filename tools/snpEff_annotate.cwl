cwlVersion: v1.0
class: CommandLineTool
id: snpEff_annotate
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/snpeff:4_3t'
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      tar -xzvf
      $(inputs.ref_tar_gz.path)
      && cwd=`pwd`
      && java -Xms2000m -Xmx8000m -jar /snpEff/snpEff.jar
      -dataDir $cwd
      -nodownload
      -t
      hg38
      $(inputs.input_vcf.path)
      | bgzip -c > $(inputs.output_basename).$(inputs.tool_name).snpEff.vcf.gz
      && tabix $(inputs.output_basename).$(inputs.tool_name).snpEff.vcf.gz
inputs:
  ref_tar_gz: { type: File, label: tar gzipped snpEff reference}
  input_vcf: { type: File,  secondaryFiles: [.tbi] }
  output_basename: string
  tool_name: string
outputs:
  output_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
  output_tbi:
    type: File
    outputBinding:
      glob: '*.vcf.gz.tbi'
