cwlVersion: v1.0
class: CommandLineTool
id: normalize_vcf
doc: "Before consensus calling, left align indel calls, break up multi-allelic calls, but leave mnps intact"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest'

baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 0
    shellQuote: true
    valueFrom: >-
      set -eo pipefail

      VCF=$(inputs.input_vcf.path)

      ${
          var cmd = " >&2 echo checking if strip flag given;";
          if (inputs.strip_info != null){
            cmd += ">&2 echo strip flag given; VCF=stripped.vcf;"
            cmd += "bcftools annotate -x " + inputs.strip_info + " " + inputs.input_vcf.path + " -o $VCF"
            cmd += " || VCF=" + inputs.input_vcf.path
          }else{
            cmd += " >&2 echo no strip flag given"
          }
          return cmd;
      }
      && bcftools norm --threads 4 -m "-any" $VCF
      | /vt/vt normalize - -n -r $(inputs.indexed_reference_fasta.path)
      | bgzip -@ 4 -c >  $(inputs.output_basename).$(inputs.tool_name).bcf_vt_norm.vcf.gz
      && tabix $(inputs.output_basename).$(inputs.tool_name).bcf_vt_norm.vcf.gz

inputs:
    input_vcf: {type: File, secondaryFiles: ['.tbi']}
    indexed_reference_fasta: {type: File, secondaryFiles: ['.fai']}
    output_basename: string
    tool_name: string
    strip_info: {type: ['null', string], doc: "If given, remove previous annotation information based on INFO file, i.e. to strip VEP info, use INFO/ANN or INFO/CSQ - check vcf"}

outputs:
  normalized_vcf:
    type: File
    outputBinding:
      glob: '*.bcf_vt_norm.vcf.gz'
    secondaryFiles: ['.tbi']
