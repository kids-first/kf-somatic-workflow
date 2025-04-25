cwlVersion: v1.2
class: CommandLineTool
id: clt_file_to_file
doc: |
  "This is a simple tool to get around a technical limitation of cavatica.
  The pickValue feature works within tools, and `outputSource`, but will not work in `outputSource`
  when used in the context of a sub-workflow. This tool gets around that limitation."
requirements:
  - class: InlineJavascriptRequirement
baseCommand: [echo, done]
inputs:
  input_file: { type: File }
outputs:
  output_file:
    type: File
    outputBinding:
      outputEval: |
        $(inputs.input_file)
