cwlVersion: v1.0
class: CommandLineTool
id: generic_rename_outputs
doc: "Simple tool to rename inputs"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 2000
    coresMin: 1

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      ${
          var out = "renamed/"
          var cmd = "mkdir " + out + ";";
          for (var i = 0; i < inputs.input_files.length; i++) {
            cmd += "cp " + inputs.input_files[i].path + " " + out + inputs.rename_to[i] + ";"
          }
          return cmd;
      }
inputs:
    input_files: {type: 'File[]', doc: "Files to rename"}
    rename_to: {type: 'string[]', doc: "List of new file names"}
outputs:
  renamed_files:
    type: 'File[]'
    outputBinding:
      glob: 'renamed/*.*'
