cwlVersion: v1.2
class: CommandLineTool
id: bcftools_reheader_samples_index
doc: |
  This tool is for renaming the sample names in a somatic BCF/VCF.GZ. That is,
  any BCF/VCF.GZ that has two samples: tumor and normal.
  
  This tool goes about the renaming in two ways:
  
  1. If BOTH old_normal_name and old_tumor_name are provided, it will construct a
  sample_list file that has the lines "old_normal_name new_normal_name" and
  "old_tumor_name new_tumor_name". This sample_list file will then be fed to
  reheader and can be used to replace the samples named old_normal_name with
  new_normal_name and old_tumor_name with new_tumor_name. This is directed,
  specific replacement. If the old names do not match any samples in the VCF,
  nothing will be changed.
  
  2. If either old_normal_name or old_tumor_name is missing, this tool will
  assume that first sample in the VCF is the normal sample and replace it with
  new_normal_name. It will assume the second sample in the VCF is the tumor
  sample and replace it with new_tumor_name. This approach is VERY DANGEROUS!
  Only use this approach if you are absolutely sure that the sample order is
  FIXED and the first sample will ALWAYS BE NORMAL!
  
  The resulting output is then indexed using BCFtools index.

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'staphb/bcftools:1.17'
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 1
  - class: InitialWorkDirRequirement
    listing:
      - entryname: sample_list.txt
        entry: |
          $(inputs.old_normal_name != null && inputs.old_tumor_name != null ? inputs.old_normal_name + " " : "")$(inputs.new_normal_name)
          $(inputs.old_tumor_name != null && inputs.old_normal_name != null ? inputs.old_tumor_name + " " : "")$(inputs.new_tumor_name)
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      bcftools reheader --samples sample_list.txt
  - position: 10
    prefix: "&&"
    shellQuote: false
    valueFrom: >-
      bcftools index
  - position: 19
    shellQuote: false
    valueFrom: >-
      $(inputs.output_filename)

inputs:
  input_vcf: { type: 'File', inputBinding: { position: 9 }, secondaryFiles: [{ pattern: '.tbi', required: false }, { pattern: '.csi', required: false }], doc: "Somatic BCF/VCF.GZ file to be reheadered with new normal and tumor sample names. Must have only two samples with the first sample being normal and the second sample being tumor" }
  output_filename: { type: 'string?', default: "reheadered.vcf.gz", inputBinding: { prefix: "--output", position: 2 }, doc: "Output name for the reheadered VCF." }

  # Somatic Samples File Reheader Options
  new_normal_name: { type: 'string', doc: "New normal sample name. Used to replace either the old_normal_name or the first sample in the input_vcf" }
  new_tumor_name: { type: 'string', doc: "New tumor sample name. Used to replace either the old_tumor_name or the second sample in the input_vcf" }
  old_tumor_name: { type: 'string?', doc: "Name of the tumor sample in the input_vcf. Will be replaced by the new_tumor_name." }
  old_normal_name: { type: 'string?', doc: "Name of the normal sample in the input_vcf. Will be replaced by the new_normal_name." }

  # BCFtools Index Options
  force: { type: 'boolean?', inputBinding: { position: 12, prefix: "--force"}, doc: "overwrite index if it already exists" }
  min_shift: { type: 'int?', inputBinding: { position: 12, prefix: "--min-shift"}, doc: "set minimal interval size for CSI indices to 2^INT [14]" }
  csi: { type: 'boolean?', inputBinding: { position: 12, prefix: "--csi"}, doc: "generate CSI-format index for VCF/BCF files [default]" }
  tbi: { type: 'boolean?', inputBinding: { position: 12, prefix: "--tbi"}, doc: "generate TBI-format index for VCF files" }
  nrecords: { type: 'boolean?', inputBinding: { position: 12, prefix: "--nrecords"}, doc: "print number of records based on existing index file" }
  stats: { type: 'boolean?', inputBinding: { position: 12, prefix: "--stats"}, doc: "print per contig stats based on existing index file" }

outputs:
  reheadered_vcf:
    type: File
    secondaryFiles: [{pattern: '.tbi', required: false}, {pattern: '.csi', required: false}]
    outputBinding:
      glob: $(inputs.output_filename)
