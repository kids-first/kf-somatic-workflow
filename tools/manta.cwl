cwlVersion: v1.0
class: CommandLineTool
id: kfdrc-manta-sv
label: Manta sv caller
doc: 'Calls structural variants.  Tool designed to pick correct run mode based on if tumor, normal, or both crams are given'
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: ${ return inputs.ram * 1000 }
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/manta:1.4.0'

baseCommand: [/manta-1.4.0.centos6_x86_64/bin/configManta.py]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      ${
        var std = " --ref " + inputs.reference.path + " --callRegions " + inputs.hg38_strelka_bed.path + " --runDir=./ && ./runWorkflow.py -m local -j " + inputs.cores + " ";
        var mv = " && mv results/variants/";
        if (typeof inputs.input_tumor_cram === 'undefined' || inputs.input_tumor_cram === null){
          var mv_cmd = mv + "diploidSV.vcf.gz " +  inputs.output_basename + ".manta.diploidSV.vcf.gz" + mv + "diploidSV.vcf.gz.tbi " + inputs.output_basename + ".manta.diploidSV.vcf.gz.tbi" + mv + "candidateSmallIndels.vcf.gz " + inputs.output_basename + ".manta.candidateSmallIndels.vcf.gz" + mv + "candidateSmallIndels.vcf.gz.tbi " + inputs.output_basename + ".manta.candidateSmallIndels.vcf.gz.tbi";
          return "--bam ".concat(inputs.input_normal_cram.path, std, mv_cmd);
        }
        else if (typeof inputs.input_normal_cram === 'undefined' || inputs.input_normal_cram === null){
          var mv_cmd = mv + "tumorSV.vcf.gz " + inputs.output_basename + ".manta.tumorSV.vcf.gz" + mv + "tumorSV.vcf.gz.tbi " + inputs.output_basename + ".manta.tumorSV.vcf.gz.tbi" + mv + "candidateSmallIndels.vcf.gz " + inputs.output_basename + ".manta.candidateSmallIndels.vcf.gz" + mv + "candidateSmallIndels.vcf.gz.tbi " + inputs.output_basename + ".manta.candidateSmallIndels.vcf.gz.tbi";
          return "--tumorBam " + inputs.input_tumor_cram.path + std + mv_cmd;
        }
        else{
          var mv_cmd = mv + "somaticSV.vcf.gz " + inputs.output_basename + ".manta.somaticSV.vcf.gz" + mv + "somaticSV.vcf.gz.tbi " + inputs.output_basename + ".manta.somaticSV.vcf.gz.tbi" + mv + "candidateSmallIndels.vcf.gz " + inputs.output_basename + ".manta.candidateSmallIndels.vcf.gz" + mv + "candidateSmallIndels.vcf.gz.tbi " + inputs.output_basename + ".manta.candidateSmallIndels.vcf.gz.tbi";
          return "--tumorBam " + inputs.input_tumor_cram.path + " --normalBam " + inputs.input_normal_cram.path + std + mv_cmd;
        }
      }

inputs:
    reference: {type: 'File', secondaryFiles: [^.dict, .fai]}
    hg38_strelka_bed: {type: 'File', secondaryFiles: [.tbi]}
    input_tumor_cram: {type: ["null", File], secondaryFiles: [.crai]}
    input_normal_cram: {type: ["null", File], secondaryFiles: [.crai]}
    cores: {type: ['null', int], default: 18}
    ram: {type: 'int?', default: 10, doc: "GB of RAM an instance must have to run the task"}
    output_basename: string
outputs:
  output_sv:
    type: File
    outputBinding:
      glob: '*SV.vcf.gz'
    secondaryFiles: [.tbi]
  small_indels:
    type: File
    outputBinding:
      glob: '*SmallIndels.vcf.gz'
    secondaryFiles: [.tbi]
