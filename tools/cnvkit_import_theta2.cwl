cwlVersion: v1.0
class: CommandLineTool
id: cnvkit_import_theta2
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4

baseCommand: [cnvkit.py, import-theta]  
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >- 
      $(inputs.tumor_cns.path)
      $(inputs.theta2_n2_results.path)
      -d ./

      mv $(inputs.output_basename).call-1.cns $(inputs.output_basename).theta2.total.cns
      
      ln -s $(inputs.output_basename).theta2.total.cns $(inputs.tumor_sample_name).cns

      cnvkit.py export seg $(inputs.tumor_sample_name).cns -o $(inputs.output_basename).theta2.total.seg

      rm $(inputs.tumor_sample_name).cns

      cnvkit.py import-theta
      $(inputs.tumor_cns.path)
      $(inputs.theta2_best_results.path)
      -d ./

      mv $(inputs.output_basename).call-1.cns $(inputs.output_basename).theta2.subclone1.cns

      ln -s $(inputs.output_basename).theta2.subclone1.cns $(inputs.tumor_sample_name).cns

      cnvkit.py export seg $(inputs.tumor_sample_name).cns -o $(inputs.output_basename).theta2.subclone1.seg

      rm $(inputs.tumor_sample_name).cns

      SC2=$(inputs.output_basename).call-2.cns;

      if [ -f "$SC2" ]; then
          mv $(inputs.output_basename).call-2.cns $(inputs.output_basename).theta2.subclone2.cns;
          ln -s $(inputs.output_basename).theta2.subclone2.cns $(inputs.tumor_sample_name).cns;
          cnvkit.py export seg $(inputs.tumor_sample_name).cns -o $(inputs.output_basename).theta2.subclone2.seg;
          rm $(inputs.tumor_sample_name).cns;
      else
        echo "second subclone file not found.  skipping!";
      fi

inputs:
  tumor_cns: File
  tumor_sample_name: string
  output_basename: string
  theta2_best_results: File
  theta2_n2_results: File

outputs:
  theta2_adjusted_cns:
    type: File
    outputBinding:
      glob: '*.theta2.total.cns'
  theta2_adjusted_seg:
    type: File
    outputBinding:
      glob: '*.theta2.total.seg'
  theta2_subclone_cns:
    type: File[]
    outputBinding:
      glob: '*.theta2.subclone*.cns'
  theta2_subclone_seg:
    type: File[]
    outputBinding:
      glob: '*.theta2.subclone*.seg'
    