cwlVersion: v1.0
class: CommandLineTool
id: gen_controlfreec_configfile
label: gen_controlfreec_configfile

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'ubuntu:18.04'
  
baseCommand: []

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      CONFIG=config.txt

      CONTENT="[general]
      
      chrLenFile = $(inputs.chr_len.path)
      
      chrFiles = ./GRCh38_everyChrs
      
      ploidy = 2
      
      maxThreads = $(inputs.threads)
      
      coefficientOfVariation = 0.062";

      printf "$CONTENT" > $CONFIG
      
      ${
        if (inputs.exome_flag == 'Y') {
          var exome = "echo \"\nwindow = 0\nreadCountThreshold = 50\nnoisyData = TRUE\nprintNA = FALSE\n\"";
          return exome;
        }
        else{
          return "echo \n";
        }
      }
      >> $CONFIG

      CONTENT="sambamba = /usr/local/bin/sambamba
      
      samtools = /usr/bin/samtools
      
      [sample]
      
      mateFile = $(inputs.tumor_bam.path)
      
      inputFormat = BAM
      
      mateOrientation = FR
      
      [control]
      
      mateFile = $(inputs.normal_bam.path)
      
      inputFormat = BAM
      
      mateOrientation = FR";

      printf "$CONTENT" >> $CONFIG
      
      ${
        if (inputs.exome_flag == 'Y') {
          var exome = "echo \"\n[target]\ncaptureRegions = " + inputs.capture_regions.path + "\n\"";
          return exome;
        }
        else{
          return "echo \n";
        }
      }
      >> $CONFIG

inputs:
  tumor_bam: File
  normal_bam: File
  chr_len: File
  threads: int
  capture_regions: {type: ['null', File], doc: "Only needed if whole exome"}
  exome_flag: string

outputs:
  config_file:
    type: File
    outputBinding:
      glob: 'config.txt'