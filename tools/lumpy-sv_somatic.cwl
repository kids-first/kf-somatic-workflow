cwlVersion: v1.0
class: CommandLineTool
id: kf-lumpy-sv-somatic
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 24000
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/lumpy:0.2.13'
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      /samtools-1.8/misc/seq_cache_populate.pl
      -root $PWD/ref_cache
      $(inputs.reference.path)
      && export REF_CACHE=$PWD/ref_cache/%2s/%2s/%s
      && export REF_PATH=$(inputs.reference.path)
      && lumpyexpress
      -B $(inputs.input_tumor_align.path),$(inputs.input_normal_align.path)
      -R $(inputs.reference.path)
      -o unsorted_results.vcf
      -v
      && cat unsorted_results.vcf | awk '$0~"^#" { print $0; next } { print $0 | "LC_ALL=C sort -k1,1V -k2,2n" }'
      | /lumpy-sv/lib/htslib/bgzip -c > $(inputs.output_basename).vcf.gz
      && /lumpy-sv/lib/htslib/tabix $(inputs.output_basename).vcf.gz

inputs:
  reference: { type: File,  secondaryFiles: [.fai], label: Fasta genome assembly with index }
  input_tumor_align:
    type: File
    secondaryFiles: |
      ${
        var path = inputs.input_tumor_align.location+".crai";
        if (inputs.input_tumor_align.nameext == ".bam"){
          path = inputs.input_tumor_align.location+".bai";
        }
        return { "location": path, "class": "File"};
      }

  input_normal_align:
      type: File
      secondaryFiles: |
        ${
          var path = inputs.input_normal_align.location+".crai";
          if (inputs.input_normal_align.nameext == ".bam"){
            path = inputs.input_normal_align.location+".bai";
          }
          return { "location": path, "class": "File"};
        }

  output_basename: string

outputs:
  output:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]
