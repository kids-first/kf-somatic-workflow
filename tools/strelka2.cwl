cwlVersion: v1.0
class: CommandLineTool
id: strelka2
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 10000
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/strelka'

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      ${
        if (inputs.extra_arg){
          var cmd = "sed  -e 's/extraVariantCallerArguments.*/extraVariantCallerArguments = "
          + inputs.extra_arg + "/' /strelka-2.9.3.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py.ini > config.ini && "
          return cmd
        }
        else{
          return "";
        }
      }
      /strelka-2.9.3.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py
      --normalBam $(inputs.input_normal_aligned.path)
      --tumorBam $(inputs.input_tumor_aligned.path)
      --ref $(inputs.reference.path)
      --callRegions $(inputs.hg38_strelka_bed.path)
      $(inputs.manta_small_indels && inputs.use_manta_small_indels ? '--indelCandidates ' + inputs.manta_small_indels.path : '')
      ${
        var arg = "--runDir=./";
        if (inputs.exome_flag == 'Y'){
          arg += " --exome"
        }
        return arg
      }
      ${
        if (inputs.extra_arg){
          return "--config config.ini"
        }
        else{
          return ""
        }
      } && ./runWorkflow.py
      -m local
      -j $(inputs.cores)


inputs:
  reference: { type: 'File', secondaryFiles: [^.dict, .fai] }
  hg38_strelka_bed: { type: 'File', secondaryFiles: [.tbi], label: gzipped bed file }
  exome_flag: { type: ['null', string], doc: "Y if exome/capture, defaults to WGS"}
  manta_small_indels: { type: 'File?', secondaryFiles: [.tbi], doc: "Small indels file from a Manta run" }
  use_manta_small_indels: { type: 'boolean?', default: false, doc: "Should the program use the small indels file? Defaults to false" }
  input_tumor_aligned:
    type: File
    secondaryFiles: |
      ${
        var dpath = self.location.replace(self.basename, "")
        if(self.nameext == '.bam'){
          return {"location": dpath+self.nameroot+".bai", "class": "File"}
        }
        else{
          return {"location": dpath+self.basename+".crai", "class": "File"}
        }
      }
    doc: "tumor BAM or CRAM"
  input_normal_aligned:
    type: File
    secondaryFiles: |
      ${
        var dpath = self.location.replace(self.basename, "")
        if(self.nameext == '.bam'){
          return {"location": dpath+self.nameroot+".bai", "class": "File"}
        }
        else{
          return {"location": dpath+self.basename+".crai", "class": "File"}
        }
      }
    doc: "normal BAM or CRAM"
  cores: {type: ['null', int], default: 16}
  extra_arg: {type: 'string?', doc: "Add special options to config file, i.e. --max-input-depth 1000"}
outputs:
  output_snv:
    type: File
    outputBinding:
      glob: 'results/variants/*.snvs.vcf.gz'
    secondaryFiles: [.tbi]
  output_indel:
    type: File
    outputBinding:
      glob: 'results/variants/*.indels.vcf.gz'
    secondaryFiles: [.tbi]

