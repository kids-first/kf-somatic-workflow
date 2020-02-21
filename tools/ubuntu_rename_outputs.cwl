cwlVersion: v1.0
class: CommandLineTool
id: ubuntu_rename_cf_outputs
doc: "Rename contrfreeec outputs"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'ubuntu:18.04'
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 2
  - class: InlineJavascriptRequirement
baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      ${
          var cmd = "";
          for (var i=0; i < inputs.input_files.length; i++){
            if (inputs.input_files[i] != null){
              var basename = inputs.input_files[i].basename;
              var fname = basename.replace("bam_", "");
              var parts = fname.split(".");
              parts.shift();
              var check = fname.substr(fname.length - 10);
              if (check == "config.txt") {
                  cmd += "cp " + inputs.input_files[i].path + " " + inputs.output_basename + ".controlfreec.config.txt;";
              } else {
              fname = inputs.output_basename + ".controlfreec." + parts.join(".");
              cmd += " cp " + inputs.input_files[i].path + " " + fname + ";";
              }
            }
            else{
              cmd += "echo A null file was detected, skipping >&2;"
            }
          }
        for (var j=0; j < inputs.input_pngs.length; j++){
            var nameroot = inputs.input_pngs[j].nameroot;
            var fname = nameroot.replace("bam_", "");
            fname = fname.replace(".txt", "");
            var parts = fname.split(".");
            parts.shift();
            fname = inputs.output_basename + ".controlfreec." + parts.join(".") + ".png";
            cmd += " cp " + inputs.input_pngs[j].path + " " + fname + ";";
            }
          return cmd;
      }

inputs:
  input_files: File[]
  input_pngs: File[]
  output_basename: string
outputs:
  ctrlfreec_baf:
    type: File?
    outputBinding:
      glob: '*.BAF.txt'
  ctrlfreec_cnv:
    type: File
    outputBinding:
      glob: '*.CNVs'
  ctrlfreec_pval:
    type: File
    outputBinding:
      glob: '*.CNVs.p.value.txt'
  ctrlfreec_config:
    type: File
    outputBinding:
      glob: '*.config.txt'
  ctrlfreec_pngs:
    type: 'File[]'
    outputBinding:
      glob: '*.png'
  ctrlfreec_bam_ratio:
    type: File
    outputBinding:
      glob: '*.ratio.txt'
  ctrlfreec_info:
    type: File
    outputBinding:
      glob: '*.info.txt'
  ctrlfreec_tumor_cpn:
    type: File
    outputBinding:
      glob: '*.sample.cpn'
  ctrlfreec_normal_cpn:
    type: File
    outputBinding:
      glob: '*.control.cpn'
