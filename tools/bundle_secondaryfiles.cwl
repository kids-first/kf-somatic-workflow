cwlVersion: v1.0
class: CommandLineTool
id: bundle_secondaryfiles 
doc: |-
  This tool takes a primary file and list of secondary files as input and passes the primary_file as
  the output with the secondary_files as secondaryFiles.
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
  - class: ShellCommandRequirement
  - class: InitialWorkDirRequirement
    listing: [$(inputs.primary_file),$(inputs.secondary_files)]
baseCommand: [echo]
inputs:
  primary_file: { type: File, doc: "Primary File" }
  secondary_files: { type: 'File[]', doc: "List of secondary files" }
outputs:
  output: 
    type: File
    outputBinding:
      glob: $(inputs.primary_file.basename)
    secondaryFiles: ${var arr = []; for (i = 0; i < inputs.secondary_files.length; i++) { if (inputs.secondary_files[i]) { arr.push(inputs.secondary_files[i].basename) } }; return arr}
