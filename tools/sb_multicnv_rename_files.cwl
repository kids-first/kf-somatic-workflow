class: CommandLineTool
cwlVersion: v1.1
id: sb_multicnv_rename_files
label: SB MultiCNV Rename Files
doc: |-
  This tool is used to rename the files in a defined way: **case_id.sample_type.caller.extension**. This information is obtained from the metadata of the file, so it needs to be defined in the file. If it is file which contains the information about the ploidy, string 'ploidy' is added to the name.

requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: ubuntu:14.04
- class: InlineJavascriptRequirement
  expressionLib:
  - |2-

    var setMetadata = function(file, metadata) {
        if (!('metadata' in file)) {
            file['metadata'] = {}
        }
        for (var key in metadata) {
            file['metadata'][key] = metadata[key];
        }
        return file
    };
    var inheritMetadata = function(o1, o2) {
        var commonMetadata = {};
        if (!o2) {
            return o1;
        };
        if (!Array.isArray(o2)) {
            o2 = [o2]
        }
        for (var i = 0; i < o2.length; i++) {
            var example = o2[i]['metadata'];
            for (var key in example) {
                if (i == 0)
                    commonMetadata[key] = example[key];
                else {
                    if (!(commonMetadata[key] == example[key])) {
                        delete commonMetadata[key]
                    }
                }
            }
            for (var key in commonMetadata) {
                if (!(key in example)) {
                    delete commonMetadata[key]
                }
            }
        }
        if (!Array.isArray(o1)) {
            o1 = setMetadata(o1, commonMetadata)
            if (o1.secondaryFiles) {
                o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)
            }
        } else {
            for (var i = 0; i < o1.length; i++) {
                o1[i] = setMetadata(o1[i], commonMetadata)
                if (o1[i].secondaryFiles) {
                    o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)
                }
            }
        }
        return o1;
    };

baseCommand: []

arguments:
  - prefix: ''
    position: 0
    valueFrom: |-
      ${ // handling no inputs
          if (inputs.input_file) {
              return 'cp'
          } else {
              return 'echo NO input given, skipping... #'
          }
      }
    shellQuote: false
  - prefix: ''
    position: 1
    valueFrom: |-
      ${
          if (inputs.input_file){
              if (inputs.input_file){
                  var case_id = inputs.input_file.metadata['case_id']
                  var tumor_id = inputs.input_file.metadata['sample_type']
                  var filename = inputs.input_file.basename
                  var dot = '.'
                  var path = inputs.input_file.basename
                  var caller = 'NotKnownCaller'
                  var ext = 'NotKnownExtension'
                  var ploidy = ''

                  // Ploidy files extensions

                  // Sclust
                  if (path.endsWith('.out_cn_summary.txt')){
                      var caller = 'sclust'
                      var ext = '.out_cn_summary.txt'
                      var ploidy = 'ploidy'
                  }
                  // Control FREEC
                  if (path.endsWith('.bam_info.txt') || path.endsWith('.info.txt')){
                      var caller = 'controlfreec'
                      var ext = '.bam_info.txt'
                      var ploidy = '.ploidy'
                  }
                  // PURPLE
                  if (path.endsWith('.purple.purity.tsv')){
                      var caller = 'purple'
                      var ext = '.purple.purity.tsv'
                      var ploidy = '.ploidy'
                  }
                  // Facets
                  if (path.endsWith('.purity_ploidy.tsv')){
                      var caller = 'purple'
                      var ext = '.purity_ploidy.tsv'
                      var ploidy = '.ploidy'
                  }


                  // Data files extensions

                  if (path.endsWith('.no_loh.cns')){
                      var caller = 'cnvkitNoLoh'
                      var ext = '.no_loh.cns'
                  }
                  if (path.endsWith('_cnv.txt')){
                      var caller = 'cnvkitNoLoh'
                      var ext = '._cnv.txt'
                  }
                  if (path.endsWith("cns")){
                      var caller = 'cnvkit'
                      var ext = '.cns'
                  }
                  //  PureCN no LOH
                  if (path.endsWith('.seg')){
                      var caller = 'purecnNoLoh'
                      var ext = '.seg'
                  }
                  // GATK called
                  if (path.endsWith('.called.seg')){
                      var caller = 'gatkCalled'
                      var ext = '.called.seg'
                  }
                  // GATK model
                  if (path.endsWith('.modelFinal.seg')){
                      var caller = 'gatkModel'
                      var ext = '.modelFinal.seg'
                  }
                  // PureCN
                  if (path.endsWith('_loh.csv')){
                      var caller = 'purecn'
                      var ext = '._loh.csv'
                  }
                  // ControlFREEC
                  if (path.endsWith('.value.txt')){
                      var caller = 'controlfreec'
                      var ext = '.value.txt'
                  }
                  // VarSimLab
                  if (path.endsWith('.CNV_results.txt')){
                      var caller = 'varsimlab'
                      var ext = '.CNV_results.txt'
                  }
                  // CNVnator
                  if (path.endsWith('calling_results.txt')){
                      var caller = 'cnvantor'
                      var ext = '.calling_results.txt'
                  }
                  //ICR96 truth
                  if (path.endsWith('.exon.csv')){
                      var caller = 'icr96'
                      var ext = '.exon.csv'
                  }
                  // Sequenza
                  if (path.endsWith('.segments.txt')){
                      var caller = 'sequenza'
                      var ext = '.segments.txt'
                  }
                  // Facets
                  if (path.endsWith('.cncf.tsv')){
                      var caller = 'facets'
                      var ext = '.cncf.tsv'
                  }
                  // SCNVSim
                  if (path.endsWith('.scnvsim.bed')){
                      var caller = 'scnvsim'
                      var ext = '.scnvsim.bed'
                  }
                  // SimulateCNVs
                  if (path.endsWith('.overlap_exon.bed')){
                      var caller = 'simulateCnvs'
                      var ext = '.overlap_exon.bed'
                  }
                  // BED
                  if (path.endsWith('.bed')){
                      var caller = 'bed'
                      var ext = '.bed'
                  }
                  // SVeto Prepare SV -> CNV
                  if (path.endsWith('.sveto-prep.cnv')){
                      var caller = 'truth'
                      var ext = '.sveto-prep.cnv'
                  }
                  // SBG Conseca CNV
                  if (path.endsWith('.conseca.tsv')){
                      var caller = 'conseca'
                      var ext = '.conseca.tsv'
                  }
                  // Sclust
                  if (path.endsWith('_allelic_states.txt')){
                      var caller = 'sclust'
                      var ext = '._allelic_states.txt'
                  }
                  // Purple Somatic
                  if (path.endsWith('.cnv.somatic.tsv')){
                      var caller = 'purpleSomatic'
                      var ext = '.cnv.somatic.tsv'
                  }
                  return  case_id + dot + tumor_id + dot + caller + ploidy + ext
              }
          }
      }
    shellQuote: false
inputs:
  input_file:
    type: File?
    inputBinding:
      position: 0
      shellQuote: false
    sbg:altPrefix: --input_file

outputs:
  output_file:
    label: output_file
    type: File?
    outputBinding:
      glob: |-
        ${
            if (inputs.input_file){
                if (inputs.input_file){
                    var case_id = inputs.input_file.metadata['case_id']
                    var tumor_id = inputs.input_file.metadata['sample_type']
                    var filename = inputs.input_file.basename
                    var dot = '.'
                    var path = inputs.input_file.basename
                    var caller = 'NotKnownCaller'
                    var ext = 'NotKnownExtension'
                    var ploidy = ''

                    // Ploidy files extensions

                    // Sclust
                    if (path.endsWith('.out_cn_summary.txt')){
                        var caller = 'sclust'
                        var ext = '.out_cn_summary.txt'
                        var ploidy = 'ploidy'
                    }
                    // Control FREEC
                    if (path.endsWith('.bam_info.txt') || path.endsWith('.info.txt')){
                        var caller = 'controlfreec'
                        var ext = '.bam_info.txt'
                        var ploidy = '.ploidy'
                    }
                    // PURPLE
                    if (path.endsWith('.purple.purity.tsv')){
                        var caller = 'purple'
                        var ext = '.purple.purity.tsv'
                        var ploidy = '.ploidy'
                    }
                    // Facets
                    if (path.endsWith('.purity_ploidy.tsv')){
                        var caller = 'purple'
                        var ext = '.purity_ploidy.tsv'
                        var ploidy = '.ploidy'
                    }


                    // Data files extensions

                    if (path.endsWith('.no_loh.cns')){
                        var caller = 'cnvkitNoLoh'
                        var ext = '.no_loh.cns'
                    }
                    if (path.endsWith('_cnv.txt')){
                        var caller = 'cnvkitNoLoh'
                        var ext = '._cnv.txt'
                    }
                    if (path.endsWith("cns")){
                        var caller = 'cnvkit'
                        var ext = '.cns'
                    }
                    //  PureCN no LOH
                    if (path.endsWith('.seg')){
                        var caller = 'purecnNoLoh'
                        var ext = '.seg'
                    }
                    // GATK called
                    if (path.endsWith('.called.seg')){
                        var caller = 'gatkCalled'
                        var ext = '.called.seg'
                    }
                    // GATK model
                    if (path.endsWith('.modelFinal.seg')){
                        var caller = 'gatkModel'
                        var ext = '.modelFinal.seg'
                    }
                    // PureCN
                    if (path.endsWith('_loh.csv')){
                        var caller = 'purecn'
                        var ext = '._loh.csv'
                    }
                    // ControlFREEC
                    if (path.endsWith('.value.txt')){
                        var caller = 'controlfreec'
                        var ext = '.value.txt'
                    }
                    // VarSimLab
                    if (path.endsWith('.CNV_results.txt')){
                        var caller = 'varsimlab'
                        var ext = '.CNV_results.txt'
                    }
                    // CNVnator
                    if (path.endsWith('calling_results.txt')){
                        var caller = 'cnvantor'
                        var ext = '.calling_results.txt'
                    }
                    //ICR96 truth
                    if (path.endsWith('.exon.csv')){
                        var caller = 'icr96'
                        var ext = '.exon.csv'
                    }
                    // Sequenza
                    if (path.endsWith('.segments.txt')){
                        var caller = 'sequenza'
                        var ext = '.segments.txt'
                    }
                    // Facets
                    if (path.endsWith('.cncf.tsv')){
                        var caller = 'facets'
                        var ext = '.cncf.tsv'
                    }
                    // SCNVSim
                    if (path.endsWith('.scnvsim.bed')){
                        var caller = 'scnvsim'
                        var ext = '.scnvsim.bed'
                    }
                    // SimulateCNVs
                    if (path.endsWith('.overlap_exon.bed')){
                        var caller = 'simulateCnvs'
                        var ext = '.overlap_exon.bed'
                    }
                    // BED
                    if (path.endsWith('.bed')){
                        var caller = 'bed'
                        var ext = '.bed'
                    }
                    // SVeto Prepare SV -> CNV
                    if (path.endsWith('.sveto-prep.cnv')){
                        var caller = 'truth'
                        var ext = '.sveto-prep.cnv'
                    }
                    // SBG Conseca CNV
                    if (path.endsWith('.conseca.tsv')){
                        var caller = 'conseca'
                        var ext = '.conseca.tsv'
                    }
                    // Sclust
                    if (path.endsWith('_allelic_states.txt')){
                        var caller = 'sclust'
                        var ext = '._allelic_states.txt'
                    }
                    // Purple Somatic
                    if (path.endsWith('.cnv.somatic.tsv')){
                        var caller = 'purpleSomatic'
                        var ext = '.cnv.somatic.tsv'
                    }
                    return  case_id + dot + tumor_id + dot + caller + ploidy + ext
                }
            }
        }
      outputEval: $(inheritMetadata(self, inputs.input_file))

$namespaces:
  sbg: https://sevenbridges.com
