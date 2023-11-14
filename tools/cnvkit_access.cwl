cwlVersion: v1.2
class: CommandLineTool
id: cnvkit_access 
doc: |
  Calculate the sequence-accessible coordinates in chromosomes from the given
  reference genome, output as a BED file.

  Many fully sequenced genomes, including the human genome, contain large regions
  of DNA that are inaccessable to sequencing. (These are mainly the centromeres,
  telomeres, and highly repetitive regions.) In the reference genome sequence
  these regions are filled in with large stretches of “N” characters. These
  regions cannot be mapped by resequencing, so CNVkit avoids them when
  calculating the antitarget bin locations.
  
  The access command computes the locations of the accessible sequence regions
  for a given reference genome based on these masked-out sequences, treating long
  spans of ‘N’ characters as the inaccessible regions and outputting the
  coordinates of the regions between them.
  
  Other known unmappable, variable, or poorly sequenced regions can be excluded
  with the -x/--exclude option. This option can be used more than once to exclude
  several BED files listing different sets of regions.
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "etal/cnvkit:0.9.3"
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)

baseCommand: [cnvkit.py, access]

inputs:
  reference_fasta: { type: File, inputBinding: { position: 9 }, doc: "Reference FASTA" }
  output_filename: { type: "string?", default: "access.bed", inputBinding: { position: 2, prefix: "--output" }, doc: "Name for output BED file." }
  min_gap_size: { type: "int?", inputBinding: { position: 2, prefix: "--min-gap-size" }, doc: "Minimum gap size between accessible sequence regions. Regions separated by less than this distance will be joined together." }
  exclude: { type: 'File?', inputBinding: { position: 2, prefix: "--exclude"}, doc: "Additional regions to exclude, in BED format. Can be used multiple times." }
  cpu: { type: 'int?', default: 8, doc: "CPUs to allocate to this task." }
  ram: { type: 'int?', default: 16, doc: "GB of RAM to allocate to this task." }

outputs:
  bed:
    type: File
    outputBinding:
      glob: $(inputs.otuput_filename)
