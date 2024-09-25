cwlVersion: v1.2
class: CommandLineTool
id: awk_chrlen_builder 
doc: >-
  This tool is an awk command that does the following:
  - Create a list of chromosome names from the input_intervals
  - Print any matching chromosome names from the reference_fai into a temp_chrs TSV
  - IF an input_variants file is provided:
    - Convert the file to a VCF
    - Create a list of chromosome names from the VCF
    - Print any matching chromosome names from the temp_chrs TSV
  - Otherwise rename the temp_chrs TSV and return it

  The reason for the optional step is that BAF calling in ControlFREEC will fail if
  a chromosome has no germline variants. Thus it is necessary to remove those chrs
  from the calling set entirely.
requirements:
  - class: DockerRequirement
    dockerPull: 'staphb/bcftools:1.20'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)
  - class: ShellCommandRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: build_chrlen.sh
        entry: |
          #!/usr/bin/env bash

          set -ex

          awk -F'\t' 'NR==FNR && $1 !~ /^@/ { out[$1]=1; next } { if (out[$1]==1) print $0 }' $(inputs.input_intervals.path) $(inputs.reference_fai.path) > temp_chrs.tsv

          if [ -s $(inputs.input_variants == null ? inputs.input_variants : inputs.input_variants.path) ]
          then
            bcftools view --no-header $(inputs.input_variants == null ? inputs.input_variants : inputs.input_variants.path) > temp_recs.vcf
            awk -F'\t' 'NR==FNR && $1 !~ /^@/ { out[$1]=1; next } { if (out[$1]==1) print $0 }' temp_recs.vcf temp_chrs.tsv > chrlen.tsv
          else
            mv temp_chrs.tsv chrlen.tsv
          fi

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: |
      /bin/bash build_chrlen.sh

inputs:
  input_intervals: { type: 'File', doc: "Intervals over which to perform calling."}
  input_variants: { type: 'File?', secondaryFiles: [{ pattern: ".tbi", required: false }], doc: "Germline Variants file" }
  reference_fai: { type: 'File', doc: "FAI index for the FASTA on which the intervals are based." }
  cpu: { type: 'int?', default: 8 }
  ram: { type: 'int?', default: 16 }

outputs:
  chrlen:
    type: File
    outputBinding:
      glob: "chrlen.tsv"
