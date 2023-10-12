cwlVersion: v1.2
class: CommandLineTool
id: clt_pass_file
doc: |
  Given a file return that file. Used to simplify some pickvalue scenarios.
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)
baseCommand: [echo, done]
inputs:
  infile: { type: 'File' }
  cpu: { type: 'int?', default: 1, doc: "CPUs to allocate to this task." }
  ram: { type: 'int?', default: 2, doc: "RAM (in GB) to allocate to this task." }
outputs:
  outfile:
    type: File
    outputBinding:
      outputEval: |
        $(inputs.infile)
