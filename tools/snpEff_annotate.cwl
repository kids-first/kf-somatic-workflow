cwlVersion: v1.0
class: CommandLineTool
id: snpEff_annotate
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
  - class: DockerRequirement
    dockerPull: 'kfdrc/snpeff:4_3t'
baseCommand: [tar, -xzvf]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.ref_tar_gz.path)
      && cwd=`pwd`
      && java -Xms2000m -Xmx8000m -jar /snpEff/snpEff.jar
      -dataDir $cwd
      -nodownload
      -t
      hg38
      $(inputs.input_vcf.path)
      | bgzip -c > $(inputs.output_basename).snpEff.vcf.gz
      && tabix $(inputs.output_basename).snpEff.vcf.gz
inputs:
  ref_tar_gz: File
  input_vcf: { type: File,  secondaryFiles: [.tbi]}
  output_basename: string
outputs:
  output:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]