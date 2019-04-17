cwlVersion: v1.0
class: CommandLineTool
id: gatk4_selectvariants
label: GATK Select PASS
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      ${
        var run_mode = inputs.mode;
        if (run_mode == 'grep' || run_mode == 'gatk'){
          var in_vcf = inputs.input_vcf.path;
          var out_vcf = inputs.output_basename + '.' + inputs.tool_name + '.PASS.vcf.gz';
          var cmd = '/gatk SelectVariants --java-options "-Xmx8000m" -V ' + in_vcf +  ' -O ' + out_vcf + ' --exclude-filtered TRUE';
          if (run_mode == 'grep'){
            cmd = 'zcat ' + in_vcf + ' | grep -E "^#|PASS" | bgzip > ' + out_vcf + '; tabix ' + out_vcf;
          }
          return cmd;
        }
        else{
          throw new Error(run_mode + ' is a not a valid mode.  Choices are gatk or grep.');
        }
      }

inputs:
  input_vcf: {type: File, secondaryFiles: [.tbi]}
  output_basename: string
  tool_name: string
  mode: {type: string, doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression"}
  
outputs:  
  pass_vcf:
    type: File
    outputBinding:
      glob: '*.PASS.vcf.gz'
    secondaryFiles: ['.tbi']

