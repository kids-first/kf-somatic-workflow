cwlVersion: v1.0
class: CommandLineTool
id: mode_selector 
doc: "Selects the appropriate input to serve as the output given the mode" 
requirements:
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/ubuntu:18.04'
  - class: InlineJavascriptRequirement
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 0
    valueFrom: >-
      set -eo pipefail

      ${
          var cmd = " >&2 echo Choosing inputs based on mode;";
          if (inputs.input_mode == 'WGS' && inputs.wgs_input == null){
            return "echo WGS run requires wgs_input >&2 && exit 1;"
          }
          else if (inputs.input_mode == 'WXS' && inputs.wxs_input == null){
            return "echo WXS run requires wxs_input >&2 && exit 1;"
          }
          return cmd;
      }
       
inputs:
  input_mode: {type: {type: enum, name: "input_mode", symbols: ["WGS", "WXS"] }, doc: "Select if this run is WGS or WXS"}
  wgs_input: {type: 'Any?', doc: "Input that should be passed when mode is WGS"}
  wxs_input: {type: 'Any?', doc: "Input that should be passed when mode is WXS"}

outputs:
  output: 
    type: Any 
    outputBinding:
      outputEval: >-
        ${
          if (inputs.input_mode == 'WGS') { return inputs.wgs_input }
          else if (inputs.input_mode == 'WXS') { return inputs.wxs_input }
        }
