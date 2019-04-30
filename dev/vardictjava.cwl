cwlVersion: v1.0
class: CommandLineTool
id: vardictjava
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 18
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'kfdrc/vardict:1.5.8'

baseCommand: [/bin/bash, -c]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail;
      export VAR_DICT_OPTS='"-Xms768m" "-Xmx16g"';
      /VarDict-1.5.8/bin/VarDict
      -G $(inputs.reference.path) -f 0.05 -th 4 --nosv -N $(inputs.output_basename)
      -b '$(inputs.input_tumor_bam.path)|$(inputs.input_normal_bam.path)'
      -z -c 1 -S 2 -E 3 -g 4 -y -F 0x700 -Q 10 -V 0.01 -Y 0 $(inputs.bed.path)
      | /VarDict-1.5.8/bin/testsomatic.R
      | /VarDict-1.5.8/bin/var2vcf_paired.pl
      -N '$(inputs.input_tumor_name)|$(inputs.input_normal_name)' -f 0.05 -M -m 4.25 > $(inputs.output_basename).result.vcf
      && cat $(inputs.output_basename).result.vcf | perl -e 'while(<>){if ($_ =~ /^#/){print $_;} else{@a = split /\t/,$_; if($a[3] =~ /[KMRYSWBVHDXkmryswbvhdx]/){$a[3] = "N";} if($a[4] =~ /[KMRYSWBVHDXkmryswbvhdx]/){$a[4] = "N";} if($a[3] ne $a[4]){print join("\t", @a);}}}'
      > $(inputs.output_basename).canonical_base_only.$(inputs.bed.nameroot).vcf
      && bgzip  $(inputs.output_basename).canonical_base_only.$(inputs.bed.nameroot).vcf
      && tabix  $(inputs.output_basename).canonical_base_only.$(inputs.bed.nameroot).vcf.gz

inputs:
  reference: {type: File, secondaryFiles: [^.dict, .fai]}
  input_tumor_bam: {type: File, secondaryFiles: [^.bai]}
  input_tumor_name: string
  input_normal_bam: {type: File, secondaryFiles: [^.bai]}
  input_normal_name: string
  output_basename: {type: string}
  bed: {type: File}
outputs:
  vardict_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]