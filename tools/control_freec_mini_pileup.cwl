cwlVersion: v1.0
class: CommandLineTool
id: controlfreec_mini_pileup
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 10000
    coresMin: $(inputs.threads)
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      ${
        if (inputs.snp_vcf == null){
          return "echo No vcf provided, skipping >&2 && exit 0;"
        }
        else{
          return "echo Creating pileup >&2;";
        }
      }

      zcat ${if (inputs.snp_vcf != null) {return inputs.snp_vcf.path}} | grep -v "#" | awk {'printf ("%s\t%s\t%s\t%s\t%s\n", $1,$2-1,$2,$4,$5)'} > snps.bed

      /opt/sambamba_0.5.9/sambamba_v0.5.9 mpileup
      -t $(inputs.threads)
      -o $(inputs.input_reads.nameroot).miniPileup $(inputs.input_reads.path)
      --samtools -f $(inputs.reference.path) -d 8000 -Q 0 -q 1 -l snps.bed

inputs:
  input_reads: {type: File, secondaryFiles: ['^.bai']}
  threads:
    type: ['null', int]
    default: 16
  reference: {type: File, secondaryFiles: [.fai]}
  snp_vcf: {type: File?, doc: "Germline vcf with sites to filter pileup on. Made optional to skip if needed during pipeline run"}
outputs:
  pileup:
    type: File?
    outputBinding:
      glob: $(inputs.input_reads.nameroot).miniPileup
