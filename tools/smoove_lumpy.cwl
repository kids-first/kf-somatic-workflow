cwlVersion: v1.0
class: CommandLineTool
id: kfdrc-smoove-lumpy-sv
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: ${return inputs.cores * 2000}
    coresMin: $(inputs.cpus)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/smoove:0.2.5'
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -euo pipefail

      mkdir TMP

      export TMPDIR=/$PWD/TMP

      ${
        if (inputs.input_tumor_align == null && inputs.input_normal_align == null){
          throw new Error('Need to provide one or both of a tumor bam/cram and normal bam/cram');
        }
        else{
          return "echo Simple alignment file input check passed";
        }
      }

      /samtools-1.9/misc/seq_cache_populate.pl
      -root $PWD/ref_cache
      $(inputs.reference.path)
      && export REF_CACHE=$PWD/ref_cache/%2s/%2s/%s
      && export REF_PATH=$(inputs.reference.path)
      && smoove
      call --name $(inputs.output_basename)
      --fasta $(inputs.reference.path)
      --processes $(inputs.cores)
      --outdir ./
      --genotype
      ${
          var args="";
          if (inputs.exclude_bed != null){
              args += " --exclude " + inputs.exclude_bed.path;
          }
          if (inputs.disable_smoove != null && inputs.disable_smoove == "Y"){
              args += " --noextrafilters";
          }
          if (inputs.duphold_flag != null && inputs.duphold_flag == "Y"){
              args += " --duphold";
          }
          if (inputs.support != null){
              args += " --support " + inputs.support;
          }
          if (inputs.input_tumor_align != null){
            args += " " + inputs.input_tumor_align.path;
          }
          if (inputs.input_normal_align != null){
            args += " " + inputs.input_normal_align.path;
          }
          return args;
      }

      tabix $(inputs.output_basename)-smoove.genotyped.vcf.gz

inputs:
  reference: { type: File,  secondaryFiles: [.fai], label: Fasta genome assembly with index }
  input_tumor_align:
    type: File?
    secondaryFiles: |
      ${
        var path = inputs.input_tumor_align.location+".crai";
        if (inputs.input_tumor_align.nameext == ".bam"){
          path = inputs.input_tumor_align.location+".bai";
        }
        return { "location": path, "class": "File"};
      }

  input_normal_align:
      type: File?
      secondaryFiles: |
        ${
          var path = inputs.input_normal_align.location+".crai";
          if (inputs.input_normal_align.nameext == ".bam"){
            path = inputs.input_normal_align.location+".bai";
          }
          return { "location": path, "class": "File"};
        }
  duphold_flag: {type: ['null', {type: enum, name: duphold_flag, symbols: ["Y", "N"] }], doc: "Run Brent P duphold and annotate DUP and DEL with depth change", default: "Y"}
  disable_smoove: {type: ['null', {type: enum, name: disable_smoove, symbols: ["Y", "N"] }], doc: "Disable smoove filtering, just do standard lumpy_filter", default: "N"}
  exclude_bed: {type: ['null', File], doc: "Bed file with regions to exlude from analysis, i.e. non-canonical chromosomes. Highly recommneded."}
  output_basename: string
  cores: {type: ['null', int], default: 16}
  support: {type: ['null', int], doc: "Min support for a variant. App default is 4."}

outputs:
  output:
    type: File
    outputBinding:
      glob: '*-smoove.genotyped.vcf.gz'
    secondaryFiles: [.tbi]
