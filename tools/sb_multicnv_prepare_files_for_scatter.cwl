class: CommandLineTool
cwlVersion: v1.1
id: sb_multicnv_prepare_files_for_scatter
label: SB MultiCNV Prepare Files For Scatter

requirements:
- class: DockerRequirement
  dockerPull: ubuntu:14.04
- class: InitialWorkDirRequirement
  listing:
  - $(inputs.input_files)
- class: InlineJavascriptRequirement


baseCommand: []

inputs:
  input_files:
    type: File[]
  override_case_id:
    type: string?
  override_sample_type:
    type: string?

outputs:
  output_files:
    type: File[]?
    outputBinding:
      outputEval: |-
        ${
            var inp=[].concat(inputs.input_files);
            var first=Array();
                for (var i=0;i<inp.length;i++){
                    first.push(inp[i]);
                    if(inp[i]) {
                        first[i].metadata=inp[i].metadata;
                        if (inputs.override_case_id != null) {
                            first[i].metadata["case_id"] = inputs.override_case_id
                        }
                        if (inputs.override_sample_type != null) {
                            first[i].metadata["sample_type"] = inputs.override_sample_type
                        }
                    }
                }
            return first;
        }

$namespaces:
  sbg: https://sevenbridges.com
