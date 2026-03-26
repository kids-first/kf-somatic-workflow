cwlVersion: v1.2
class: CommandLineTool
id: samtools-mini-pileup
doc: "Create mini pileup for ControlFreeC"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.20-multi-arch'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: $(inputs.threads)

baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail;

      $(inputs.snp_vcf == null ? "echo No vcf provided, skipping >&2 && exit 0" : "echo Creating pileup >&2");

      pigz -dc $(inputs.snp_vcf != null ? inputs.snp_vcf.path : "") | grep -v "#" | awk {'printf ("%s\t%s\t%s\t%s\t%s\n", $1,$2-1,$2,$4,$5)'} > snps.bed;
  - position: 10
    shellQuote: false
    valueFrom: >-
      split --number=l/$(inputs.threads) snps.bed --additional-suffix=".bed" -d regions && ls regions*.bed > bed_list.txt;

      cat bed_list.txt
      | xargs -IFN -P $(inputs.threads)
      samtools mpileup
      -d 8000
      -Q 0
      -q 1
      -l FN
      -o FN.miniPileup
  - position: 20
    shellQuote: false
    valueFrom: >-
      && cat xa*.miniPileup > $(inputs.input_reads.nameroot).miniPileup

inputs:
  input_reads: {type: File, secondaryFiles: [{pattern: ".bai", required: false}, {pattern: "^.bai", required: false}, {pattern: ".crai", required: false},
      {pattern: "^.crai", required: false}],
      inputBinding: { position: 19 } }
  threads:
    type: ['null', int]
    default: 8
  reference: {type: File, secondaryFiles: [.fai],
    inputBinding: { position: 11, prefix: "-f"}}
  snp_vcf: {type: 'File?', doc: "Germline vcf with sites to filter pileup on. Made optional to skip if needed during pipeline run"}
outputs:
  pileup:
    type: File?
    outputBinding:
      glob: $(inputs.input_reads.nameroot).miniPileup
