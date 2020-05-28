cwlVersion: v1.0
class: CommandLineTool
id: mode_defaults
doc: "Selects the appropriate defaults based on given the mode"
requirements:
  - class: InlineJavascriptRequirement
baseCommand: echo
arguments:
  - position: 0
    valueFrom: >-
      Selecting $(inputs.input_mode) default values

inputs:
  input_mode: {type: {type: enum, name: "input_mode", symbols: ["WGS", "WXS"] }, doc: "Select if this run is WGS or WXS"}
  exome_flag: {type: 'string?', doc: "Whether to run in exome mode for callers."}
  cnvkit_wgs_mode: {type: 'string?', doc: "Entering Y will run cnvkit in WGS mode, otherwise it will run in hybrid mode. Defaults to Y in wgs mode."}
  i_flag: {type: 'string?', doc: "Flag to intersect germline calls on padded regions.  Use N if you want to skip this. Defaults to N in WGS mode."}
  lancet_padding: {type: 'int?', doc: "Recommend 0 if interval file padded already, half window size if not. Recommended: 0 for WXS; 300 for WGS"}
  lancet_window: {type: 'int?', doc: "window size for lancet.  default is 600, recommend 500 for WGS, 600 for exome+"}
  vardict_padding: {type: 'int?', doc: "Padding to add to input intervals, recommend 0 if intervals already padded such as in WXS, 150 if not such as in WGS"}

outputs:
  out_exome_flag:
    type: 'string?'
    outputBinding:
      outputEval: >-
        ${
          if (inputs.exome_flag) { return inputs.exome_flag }
          else if (inputs.input_mode == 'WGS') { return 'N' }
          else if (inputs.input_mode == 'WXS') { return 'Y' }
        }
  out_cnvkit_wgs_mode:
    type: 'string?'
    outputBinding:
      outputEval: >-
        ${
          if (inputs.cnvkit_wgs_mode) { return inputs.cnvkit_wgs_mode }
          else if (inputs.input_mode == 'WGS') { return 'Y' }
          else if (inputs.input_mode == 'WXS') { return 'N' }
        }
  out_i_flag:
    type: 'string?'
    outputBinding:
      outputEval: >-
        ${
          if (inputs.i_flag) { return inputs.i_flag }
          else if (inputs.input_mode == 'WGS') { return 'N' }
          else if (inputs.input_mode == 'WXS') { return null }
        }
  out_lancet_padding:
    type: 'int?'
    outputBinding:
      outputEval: >-
        ${
          if (inputs.lancet_padding) { return inputs.lancet_padding }
          else if (inputs.input_mode == 'WGS') { return 300 }
          else if (inputs.input_mode == 'WXS') { return 0 }
        }
  out_lancet_window:
    type: 'int?'
    outputBinding:
      outputEval: >-
        ${
          if (inputs.lancet_window) { return inputs.lancet_window }
          else if (inputs.input_mode == 'WGS') { return 600 }
          else if (inputs.input_mode == 'WXS') { return 600 }
        }
  out_vardict_padding:
    type: 'int?'
    outputBinding:
      outputEval: >-
        ${
          if (inputs.vardict_padding) { return inputs.vardict_padding }
          else if (inputs.input_mode == 'WGS') { return 150 }
          else if (inputs.input_mode == 'WXS') { return 0 }
        }
