cwlVersion: v1.0
class: CommandLineTool
id: vardictjava
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: ${return inputs.ram * 1000}
    coresMin: $(inputs.cpus)
  - class: DockerRequirement
    dockerPull: 'kfdrc/vardict:1.7.0'

baseCommand: [/bin/bash, -c]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail;
      ${
        var ram = Math.floor(inputs.ram/1.074 - 1);
        var exp_cmd = "export VAR_DICT_OPTS='\"-Xms768m\" \"-Xmx" + ram + "g\"';";
        return exp_cmd;
      }
      /VarDict-1.7.0/bin/VarDict
      -G $(inputs.reference.path) -f $(inputs.min_vaf) -th $(inputs.cpus) --nosv -N $(inputs.output_basename)
      -b '$(inputs.input_tumor_bam.path)|$(inputs.input_normal_bam.path)'
      -z -c 1 -S 2 -E 3 -g 4 -F 0x700 -Q 10 -V 0.01 -x $(inputs.padding) $(inputs.bed.path) > vardict_results.txt
      && cat vardict_results.txt | /VarDict-1.7.0/bin/testsomatic.R > vardict_r_test_results.txt
      && cat vardict_r_test_results.txt | /VarDict-1.7.0/bin/var2vcf_paired.pl
      -N '$(inputs.input_tumor_name)|$(inputs.input_normal_name)' -f $(inputs.min_vaf) -M -m 4.25 > $(inputs.output_basename).result.vcf
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
  cpus: {type: ['null', int], default: 9}
  ram: {type: ['null', int], default: 18, doc: "In GB"}
  padding: {type: ['null', int], doc: "Padding to add to input intervals, recommened 0 if intervals already padded, 150 if not", default: 150}
  min_vaf:
    type: ['null', float]
    doc: "Recommend 0.05"
    default: 0.05
  output_basename: {type: string}
  bed: {type: File}
outputs:
  vardict_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]