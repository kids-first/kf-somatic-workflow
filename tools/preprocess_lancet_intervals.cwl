cwlVersion: v1.0
class: CommandLineTool
id: bcftools_reheader_vcf
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/bedops:2.4.36'
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: InlineJavascriptRequirement
baseCommand: [/bin/bash, -c]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      ${
          var flag = 0
          var cmd = "";
          var bed = []
          if(inputs.strelka2_vcf != null){
              flag = 1;
              cmd += "gunzip -c " + inputs.strelka2_vcf.path + " > " + inputs.strelka2_vcf.basename + ";";
              cmd += "vcf2bed --insertions < " + inputs.strelka2_vcf.basename + " | cut -f 1-3 > strelka2.insertions.bed;";
              bed.push("strelka2.insertions.bed");
              cmd += "vcf2bed --deletions < " + inputs.strelka2_vcf.basename + " | cut -f 1-3 > strelka2.deletions.bed;";
              bed.push("strelka2.deletions.bed");
              cmd += "vcf2bed --snvs < " + inputs.strelka2_vcf.basename + " | cut -f 1-3 > strelka2.snvs.bed;";
              bed.push("strelka2.snvs.bed");
          }
          if(inputs.mutect2_vcf != null){
              flag = 1;
              cmd += "gunzip -c " + inputs.mutect2_vcf.path + " > " + inputs.mutect2_vcf.basename + ";";
              cmd += "vcf2bed --insertions < " + inputs.mutect2_vcf.basename +  " | cut -f 1-3 > mutect2.insertions.bed;";
              bed.push("mutect2.insertions.bed");
              cmd += "vcf2bed --deletions < " + inputs.mutect2_vcf.basename + " | cut -f 1-3 > mutect2.deletions.bed;";
              bed.push("mutect2.deletions.bed");
              cmd += "vcf2bed --snvs < " + inputs.mutect2_vcf.basename +  " | cut -f 1-3 > mutect2.snvs.bed;";
              bed.push("mutect2.snvs.bed");
          }
        if(flag == 0){
            cmd += "echo \"No input vcfs found to convert.  Returning ref bed\"; cp " + inputs.ref_bed.path + " " + inputs.output_basename + ".lancet_intvervals.bed;";
        }
        else{
            cmd += "cat " + bed.join(" ") + " " + inputs.ref_bed.path + " | bedtools sort | bedtools merge > " + inputs.output_basename + ".lancet_intvervals.bed;";
        }
        return cmd;
      }

inputs:
  strelka2_vcf: {type: ['null', File], doc: "PASS vcf from strelka2 run for the sample to be analyzed"}
  mutect2_vcf: {type: ['null', File], doc: "PASS vcf from mutect2 run for the sample to be analyzed"}
  ref_bed: {type: File, doc: "Exonic bed file recommended.  Will be augmented with strelka2 and/or mutect2 sites if supplied"}
  output_basename: string
outputs:
  run_bed:
    type: File
    outputBinding:
      glob: '*.lancet_intvervals.bed'
