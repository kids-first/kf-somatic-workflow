cwlVersion: v1.0
class: ExpressionTool
id: expression_preparerg
doc: |
  "This is a simple tool to get around a technical limitation of cavatica.
  The pickValue feature works within tools, and `outputSource`, but will not work in `outputSource`
  when used in the context of a sub-workflow. This tool gets around that limitation."
requirements:
  - class: InlineJavascriptRequirement
inputs:
  input_file: { type: File }
outputs:
  output: { type: File }

expression: |
  ${ return {output: inputs.input_file } }