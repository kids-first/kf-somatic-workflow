cwlVersion: v1.2
class: ExpressionTool
id: aa_run_if_wgs
requirements:
  - class: InlineJavascriptRequirement

inputs:
  mosek_license_file: 'File?'
  wgs_or_wxs: {type: {type: enum, name: wgs_or_wxs, symbols: ["WGS", "WXS"]}, doc: "Select\
      \ if this run is WGS or WXS"}

outputs:
  aa_mosek_license_file:
    type: 'File?'

expression: |
  ${
    if (inputs.wgs_or_wxs == 'WGS' && inputs.mosek_license_file){
        return {'aa_mosek_license_file': inputs.mosek_license_file };
    }
    else {
      return {'aa_mosek_license_file': null};
    }
  }
