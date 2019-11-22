cwlVersion: v1.0
class: CommandLineTool
id: kfdrc-vep-somatic-annotate-maf
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 24000
    coresMin: 16
  - class: DockerRequirement
    dockerPull: 'kfdrc/vep:r93_v2'
baseCommand: [tar, -xzf ]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.cache.path)
      && gunzip -c $(inputs.input_vcf.path) > input_file.vcf
      && perl /vcf2maf/vcf2maf.pl
      --input-vcf input_file.vcf
      --output-maf $(inputs.output_basename).$(inputs.tool_name).vep.maf
      --filter-vcf 0
      --vep-path /ensembl-vep/
      --vep-data $PWD
      --vep-forks 16
      --ncbi-build $(inputs.ref_build)
      --cache-version $(inputs.cache_version)
      --ref-fasta $(inputs.reference.path)
      --tumor-id $(inputs.tumor_id)
      --normal-id $(inputs.normal_id)
      && mv input_file.vep.vcf $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf
      && /ensembl-vep/htslib/bgzip $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf
      && /ensembl-vep/htslib/tabix $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf.gz

inputs:
  reference: {type: File,  secondaryFiles: [.fai], label: Fasta genome assembly with index}
  input_vcf:
    type: File
    secondaryFiles: [.tbi]
  output_basename: string
  tumor_id: string
  normal_id: string
  tool_name: string
  cache: { type: File, label: tar gzipped cache from ensembl/local converted cache }
  cache_version: {type: ['null', int], doc: "Version being used, should match build version", default: 93}
  ref_build: {type: ['null', string], doc: "Genome ref build used, should line up with cache.", default: "GRCh38" }


outputs:
  output_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
  output_tbi:
    type: File  
    outputBinding:
      glob: '*.vcf.gz.tbi'
  output_maf:
    type: File
    outputBinding:
      glob: '*.maf'
  warn_txt:
    type: ["null", File]
    outputBinding:
      glob: '*.txt'
