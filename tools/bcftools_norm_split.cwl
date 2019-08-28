cwlVersion: v1.0
class: CommandLineTool
id: bcftools_norm_split_vcf
doc: "Normalized a vcf by splitting multi-allelic locations, then splits into snps+mnps and indels files"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/bvcftools:latest'
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 1
  - class: InlineJavascriptRequirement
baseCommand: [bcftools,norm,-c,w,-m,-any,-Oz,-f]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.reference.path) $(inputs.input_vcf.path) >  $(inputs.input_vcf.nameroot.replace('.vcf','')).bcfNorm.vcf.gz

      ${
          var snp_cmd = "bcftools view --exclude-types ref " + inputs.input_vcf.nameroot.replace('.vcf','') + ".bcfNorm.vcf.gz";
          var indel_cmd = "bcftools view --exclude-types ref,bnd " + inputs.input_vcf.nameroot.replace('.vcf','') + ".bcfNorm.vcf.gz";
          if (inputs.grep_v_phrase != null){
              snp_cmd += " | grep -E \"^#|" + inputs.grep_v_phrase + "\"";
              indel_cmd += " | grep -E \"^#|" + inputs.grep_v_phrase + "\"";
          }
          snp_cmd += " | bcftools view --types snps,mnps -Oz - > " + inputs.input_vcf.nameroot.replace('.vcf','') + ".bcfNorm.snv_mnv.vcf.gz; ";
          indel_cmd += " | bcftools view --types indels -Oz - > " + inputs.input_vcf.nameroot.replace('.vcf','') + ".bcfNorm.indels.vcf.gz;";
          return snp_cmd + indel_cmd;
      }

      ls *.vcf.gz | xargs -IFN tabix FN

inputs:
  input_vcf: {type: File, secondaryFiles: ['.tbi']}
  reference: {type: File, secondaryFiles: ['.fai']}
  grep_v_phrase: {type: ['null', string], doc: "Specific to vardict, to filter on somatic calls"}
outputs:
  normalized_vcf:
    type: File
    outputBinding:
      glob: '*.bcfNorm.vcf.gz'
    secondaryFiles: [.tbi]
  snv_mnv_vcf:
    type: File
    outputBinding:
      glob: '*.bcfNorm.snv_mnv.vcf.gz'
    secondaryFiles: [.tbi]
  indelvcf:
    type: File
    outputBinding:
      glob: '*.bcfNorm.indels.vcf.gz'
    secondaryFiles: [.tbi]
