cwlVersion: v1.2
class: CommandLineTool
id: runtime_validator
doc: "Additional validation for user inputs."
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      echo "Validating Inputs"

inputs:
  is_wgs: {type: 'boolean', doc: "True for WGS; False for WXS."}
  vardict: {type: 'boolean', doc: "Set to false to disable Vardict. Warning: Vardict is required to run Theta2!"}
  mutect2: {type: 'boolean', doc: "Set to false to disable Mutect2. Warning: Mutect2 or Strelka2 is required to run Lancet in WGS mode!"}
  strelka2: {type: 'boolean', doc: "Set to false to disable Strelka2. Warning: Strelka2 or Mutect2 is required to run Lancet in WGS mode!"}
  lancet: {type: 'boolean', doc: "Set to false to disable Lancet."}
  controlfreec: {type: 'boolean', doc: "Set to false to disable ControlFreeC."}
  cnvkit: {type: 'boolean', doc: "Set to false to disable CNVkit. Warning: CNVkit is required to run both Amplicon Architect and Theta2!"}
  amplicon_architect: {type: 'boolean', doc: "Set to false to disable Amplicon Architect. Always disabled for WXS!"}
  theta2: {type: 'boolean', doc: "Set to false to disable Theta2."}
  manta: {type: 'boolean', doc: "Set to false to disable Manta."}
  gatk_cnv: {type: 'boolean', doc: "Set to false to disable GATK CNV."}
  mosek_present: {type: 'boolean', doc: "Mosek file is present. Warning: Mosek file is required to run Amplicon Architect!"}
  pon_present: {type: 'boolean', doc: "GATK CNV PON file is present. Warning: Panel of Normals is required to run GATK CNV!"}
  exome_flag: {type: 'string?', doc: "Whether to run in exome mode for callers."}
  cnvkit_wgs_mode: {type: 'string?', doc: "Entering Y will run cnvkit in WGS mode, otherwise it will run in hybrid mode. Defaults to Y in wgs mode."}
  i_flag: {type: 'string?', doc: "Flag to intersect germline calls on padded regions.  Use N if you want to skip this. Defaults to N in WGS mode."}
  lancet_padding: {type: 'int?', doc: "Recommend 0 if interval file padded already, half window size if not. Recommended: 0 for WXS; 300 for WGS"}
  lancet_window: {type: 'int?', doc: "window size for lancet.  default is 600, recommend 500 for WGS, 600 for exome+"}
  vardict_padding: {type: 'int?', doc: "Padding to add to input intervals, recommend 0 if intervals already padded such as in WXS, 150 if not such as in WGS"}

outputs:
  out_wgs:
    type: 'boolean'
    outputBinding:
      outputEval: >-
        $(inputs.is_wgs)
  run_vardict:
    type: 'boolean'
    outputBinding:
      outputEval: >-
        $(inputs.vardict)
  run_mutect2:
    type: 'boolean'
    outputBinding:
      outputEval: >-
        $(inputs.mutect2)
  run_strelka2:
    type: 'boolean'
    outputBinding:
      outputEval: >-
        $(inputs.strelka2)
  run_lancet:
    type: 'boolean'
    outputBinding:
      outputEval: >-
        ${
          if (!inputs.lancet || !inputs.is_wgs) {
            return inputs.lancet;
          }
          if (!(inputs.mutect2 || inputs.strelka2)) {
            throw new Error('Running Lancet in WGS Mode requires either Mutect2 or Strelka2!');
          }
          return true;
        }
  run_controlfreec:
    type: 'boolean'
    outputBinding:
      outputEval: >-
        $(inputs.controlfreec)
  run_cnvkit:
    type: 'boolean'
    outputBinding:
      outputEval: >-
        $(inputs.cnvkit)
  run_amplicon_architect:
    type: 'boolean'
    outputBinding:
      outputEval: >-
        ${
          if (!inputs.amplicon_architect || !inputs.is_wgs) {
            return false;
          }
          if (!inputs.cnvkit) {
            throw new Error('Amplicon Architect cannot be run without CNVkit!');
          }
          if (!inputs.mosek_present) {
            throw new Error('Amplicon Architect cannot be run without Mosek file!');
          }
          return true;
        }
  run_theta2:
    type: 'boolean'
    outputBinding:
      outputEval: >-
        ${
          if (!inputs.theta2) {
            return false;
          }
          if (!(inputs.vardict && inputs.cnvkit)) {
            throw new Error('Running Theta2 requires both Vardict and CNVkit!');
          }
          return true;
          }
  run_manta:
    type: 'boolean'
    outputBinding:
      outputEval: >-
        $(inputs.manta)
  run_gatk_cnv:
    type: 'boolean'
    outputBinding:
      outputEval: >-
        ${
          if (!inputs.gatk_cnv) {
            return false;
          }
          if (!inputs.pon_present) {
            throw new Error('Running GATK CNV requires a Panel of Normals!');
          }
          return true;
        }
  out_exome_flag:
    type: 'string?'
    outputBinding:
      outputEval: >-
        ${
          if (inputs.exome_flag) { return inputs.exome_flag }
          else if (inputs.is_wgs) { return 'N' }
          else if (!inputs.is_wgs) { return 'Y' }
        }
  out_cnvkit_wgs_mode:
    type: 'string?'
    outputBinding:
      outputEval: >-
        ${
          if (inputs.cnvkit_wgs_mode) { return inputs.cnvkit_wgs_mode }
          else if (inputs.is_wgs) { return 'Y' }
          else if (!inputs.is_wgs) { return 'N' }
        }
  out_i_flag:
    type: 'string?'
    outputBinding:
      outputEval: >-
        ${
          if (inputs.i_flag) { return inputs.i_flag }
          else if (inputs.is_wgs) { return 'N' }
          else if (!inputs.is_wgs) { return null }
        }
  out_lancet_padding:
    type: 'int?'
    outputBinding:
      outputEval: >-
        ${
          if (inputs.lancet_padding) { return inputs.lancet_padding }
          else if (inputs.is_wgs) { return 300 }
          else if (!inputs.is_wgs) { return 0 }
        }
  out_lancet_window:
    type: 'int?'
    outputBinding:
      outputEval: >-
        ${
          if (inputs.lancet_window) { return inputs.lancet_window }
          else if (inputs.is_wgs) { return 500 }
          else if (!inputs.is_wgs) { return 600 }
        }
  out_vardict_padding:
    type: 'int?'
    outputBinding:
      outputEval: >-
        ${
          if (inputs.vardict_padding) { return inputs.vardict_padding }
          else if (inputs.is_wgs) { return 150 }
          else if (!inputs.is_wgs) { return 0 }
        }
