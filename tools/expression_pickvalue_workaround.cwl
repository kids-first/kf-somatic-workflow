cwlVersion: v1.0
class: ExpressionTool
id: expression_preparerg
requirements:
  - class: InlineJavascriptRequirement
inputs:
  input_file: { type: File }
outputs:
  output: { type: File }

expression: |
  ${ return {output: inputs.input_file } }