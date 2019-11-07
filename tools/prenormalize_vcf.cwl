cwlVersion: v1.0
class: CommandLineTool
id: prenormalize_vcf
doc: "Before consensus calling, left align indel calls, break up multi-allelic calls, but leave mnps intact"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'migbro/vcfutils:latest'

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      VCF=$(inputs.input_vcf.path)

      ${
          var cmd = "echo checking if strip flag given;";
          if (inputs.strip_info != null){
            cmd += "/vt/vt rminfo -t " + inputs.strip_info + " " + inputs.input_vcf.path + " -o stripped.vcf;"
            cmd += "VCF=stripped.vcf;echo strip flag given;" 
          }else{
            cmd += "echo no strip flag given;"
          }
          return cmd;
      }

      bcftools norm -m '-any' $VCF > $(inputs.output_basename).$(inputs.tool_name).bcf_norm.vcf
      && /vt/vt normalize -r $(inputs.indexed_reference_fasta.path) $(inputs.output_basename).$(inputs.tool_name).bcf_norm.vcf
      > $(inputs.output_basename).$(inputs.tool_name).bcf_vt_norm.vcf
      && bgzip $(inputs.output_basename).$(inputs.tool_name).bcf_vt_norm.vcf
      && tabix $(inputs.output_basename).$(inputs.tool_name).bcf_vt_norm.vcf.gz

inputs:
    input_vcf: {type: File, secondaryFiles: ['.tbi']}
    indexed_reference_fasta: {type: File, secondaryFiles: ['.fai']}
    output_basename: string
    tool_name: string
    strip_info: {type: ['null', string], doc: "If given, remove previous annotation information based on INFO file, i.e. to strip VEP info, use ANN"}

outputs:
  normalized_vcf:
    type: File
    outputBinding:
      glob: '*.bcf_vt_norm.vcf.gz'
    secondaryFiles: ['.tbi']
