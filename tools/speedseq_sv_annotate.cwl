cwlVersion: v1.0
class: CommandLineTool
id: kfdrc-speedseq-sv-annotate
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 24000
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/speedseq:latest'
baseCommand: [/speedseq/src/samtools-1.8/misc/seq_cache_populate.pl ]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -root $PWD/ref_cache
      $(inputs.reference.path)
      && export REF_CACHE=$PWD/ref_cache/%2s/%2s/%s
      && export REF_PATH=$(inputs.reference.path)
      && /speedseq/bin/speedseq sv
      -B $(inputs.input_align.path)
      -R $(inputs.reference.path)
      -t 8
      -o $(inputs.output_basename)
      -v

inputs:
  reference: { type: File,  secondaryFiles: [.fai], label: Fasta genome assembly with index }
  input_align:
    type: File
    secondaryFiles: |
      ${
        var path = inputs.input_align.location+".crai";
        if (inputs.input_align.nameext == ".bam"){
          path = inputs.input_align.location+".bai";
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
