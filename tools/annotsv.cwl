class: CommandLineTool
cwlVersion: v1.2
id: annotsv
requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: pgc-images.sbgenomics.com/d3b-bixu/annotsv:3.1.1
- class: InlineJavascriptRequirement
  expressionLib:
  - |2-

    var setMetadata = function(file, metadata) {
        if (!('metadata' in file))
            file['metadata'] = metadata;
        else {
            for (var key in metadata) {
                file['metadata'][key] = metadata[key];
            }
        }
        return file
    };

    var inheritMetadata = function(o1, o2) {
        var commonMetadata = {};
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
        }
        if (!Array.isArray(o1)) {
            o1 = setMetadata(o1, commonMetadata)
        } else {
            for (var i = 0; i < o1.length; i++) {
                o1[i] = setMetadata(o1[i], commonMetadata)
            }
        }
        return o1;
    };

baseCommand: []
arguments:
- position: 0
  prefix: ''
  shellQuote: false
  valueFrom: |
    tar -xf $(inputs.annotations_dir_tgz.path) --use-compress-program=pigz && AnnotSV -annotationsDir ./AnnotSV
- position: 99
  prefix: ''
  shellQuote: false
  valueFrom: |
    1>&2

inputs:
  annotations_dir_tgz: { type: 'File', doc: "TAR.GZ'd Directory containing annotations" }
  annotation_mode:
    type:
      - 'null'
      - type: enum
        name: annotaiton_mode
        symbols: ["both","full","split"]
    inputBinding:
      prefix: "-annotationMode"
      position: 1
    doc: "Description of the types of lines produced by AnnotSV"
  benign_af: { type: 'float?', inputBinding: { prefix: "-benignAF", position: 1 }, doc: "Allele frequency threshold to select the benign SV in the data sources" }
  candidate_genes_file: { type: 'File?', inputBinding: { prefix: "-candidateGenesFile", position: 1 }, doc: "Path of a file containing the candidate genes of the user (gene names can be space-separated, tabulation-separated, or line-break-separated)",  "sbg:fileTypes": "BED, CSV, TSV, TXT"}
  candidate_genes_filtering:
    type:
      - 'null'
      - type: enum
        name: candidate_genes_filtering
        symbols: ["0","1"]
    inputBinding:
      prefix: "-candidateGenesFiltering"
      position: 1
    doc: "To select only the SV annotations ('split' and 'full') overlapping a gene from the candidateGenesFile"
  candidate_snv_indel_files: { type: 'File[]?', secondaryFiles: [{pattern: ".tbi", required: false}], inputBinding: { prefix: "-candidateSnvIndelFiles", position: 1 }, doc: "Path of the filtered VCF input file(s) with SNV/indel coordinates for compound heterozygotes report (optional)", "sbg:fileTypes": "VCF, VCF.GZ" }
  candidate_snv_indel_samples: { type: 'string[]?', inputBinding: { prefix: "-candidateSnvIndelSamples", position: 1 }, doc: "To specifiy the sample names from the VCF files defined from the -candidateSnvIndelFiles option" }
  genome_build:
    type:
      - 'null'
      - type: enum
        name: genome_build
        symbols: ["GRCh37","GRCh38","mm9","mm10"]
    inputBinding:
      prefix: "-genomeBuild"
      position: 1
    doc: "Genome build used"
  hpo: { type: 'string[]?', inputBinding: { prefix: "-hpo", position: 1 }, doc: "HPO term(s) describing the phenotype of the individual being investigated. (e.g. HP:0001156)" }
  include_ci:
    type:
      - 'null'
      - type: enum
        name: include_ci
        symbols: ["0","1"]
    inputBinding:
      prefix: "-includeCI"
      position: 1
    doc: "To expand the 'start' and 'end' SV positions with the VCF confidence intervals (CIPOS, CIEND) around the breakpoints"
  metrics:
    type:
      - 'null'
      - type: enum
        name: metrics
        symbols: ["us","fr"]
    inputBinding:
      prefix: "-metrics"
      position: 1
    doc: "Changing numerical values from frequencies to us or fr metrics (e.g. 0.2 or 0,2)"
  min_total_number: { type: 'int?', inputBinding: { prefix: "-minTotalNumber", position: 1 }, doc: "Minimum number of individuals tested to consider a benign SV for the ranking. Range values: [100-1000]" }
  output_dir: { type: 'string?', default: '.', inputBinding: { prefix: "-outputDir", position: 1 }, doc: "String representing path to directory where outputs will be placed" }
  output_file: { type: 'string?', inputBinding: { prefix: "-outputFile", position: 1 }, doc: "Name for output file" }
  overlap: { type: 'int?', inputBinding: { prefix: "-overlap", position: 1 }, doc: "Minimum overlap (%) between user features (User BED) and the annotated SV to be reported. Range values: [0-100]" }
  overwrite:
    type:
      - 'null'
      - type: enum
        name: overwrite
        symbols: ["0","1"]
    inputBinding:
      prefix: "-overwrite"
      position: 1
    doc: "To overwrite existing output results"
  promoter_size: { type: 'int?', inputBinding: { prefix: "-promoterSize", position: 1 }, doc: "Number of bases upstream from the transcription start site" }
  rank_filtering: { type: 'string?', inputBinding: { prefix: "-rankFiltering", position: 1 }, doc: "To select the SV of a user-defined specific class (from 1 to 5). Values: use comma separated class values, or use a dash to denote a range of values (e.g.: '3,4,5' or '3-5')." }
  reciprocal:
    type:
      - 'null'
      - type: enum
        name: reciprocal
        symbols: ["0","1"]
    inputBinding:
      prefix: "-reciprocal"
      position: 1
    doc: "Use of a reciprocal overlap between SV and user features (only for annotations with features overlapping the SV)"
  re_report:
    type:
      - 'null'
      - type: enum
        name: re_report
        symbols: ["0","1"]
    inputBinding:
      prefix: "-REreport"
      position: 1
    doc: "Create a report to link the annotated SV and the overlapped regulatory elements (coordinates and sources)"
  re_select_1:
    type:
      - 'null'
      - type: enum
        name: re_select_1
        symbols: ["0","1"]
    inputBinding:
      prefix: "-REselect1"
      position: 1
    doc: "To report only the morbid, HI, TS, candidate and phenotype matched genes"
  re_select_2:
    type:
      - 'null'
      - type: enum
        name: re_select_2
        symbols: ["0","1"]
    inputBinding:
      prefix: "-REselect2"
      position: 1
    doc: "To report only the genes not present in 'Gene_name'"
  samples_id_bed_col: { type: 'int?', inputBinding: { prefix: "-samplesidBEDcol", position: 1 }, doc: "Number of the column reporting the samples ID for which the SV was called (if the input SV file is a BED)" }
  snv_indel_files: { type: 'File[]?', secondaryFiles: [{pattern: ".tbi", required: false}], inputBinding: { prefix: "-snvIndelFiles", position: 1 }, doc: "VCF input file(s) with SNV/indel coordinates used for false positive discovery. Use counts of the homozygous and heterozygous variants", "sbg:fileTypes": "VCF, VCF.GZ" }
  snv_indel_pass:
    type:
      - 'null'
      - type: enum
        name: snv_indel_pass
        symbols: ["0","1"]
    inputBinding:
      prefix: "-snvIndelPASS"
      position: 1
    doc: "To only use variants from VCF input files that passed all filters during the calling (FILTER column value equal to PASS)"
  snv_indel_samples: { type: 'string[]?', inputBinding: { prefix: "-snvIndelSamples", position: 1 }, doc: "To specify the sample names from the VCF files defined from the -snvIndelFiles option" }
  sv_input_file: { type: 'File', secondaryFiles: [{pattern: ".tbi", required: false}], inputBinding: { prefix: "-SVinputFile", position: 1 }, doc: "Path of the input file (VCF or BED) with SV coordinates. Gzipped VCF file is supported", "sbg:fileTypes": "BED, VCF, VCF.GZ" }
  sv_input_info:
    type:
      - 'null'
      - type: enum
        name: sv_input_info
        symbols: ["0","1"]
    inputBinding:
      prefix: "-SVinputInfo"
      position: 1
    doc: "To extract the additional SV input fields and insert the data in the outputfile"
  sv_min_size: { type: 'int?', inputBinding: { prefix: "-SVminSize", position: 1 }, doc: "SV minimum size (in bp)" }
  svt_bed_col: { type: 'int?', inputBinding: { prefix: "-svtBEDcol", position: 1 }, doc: "Number of the column describing the SV type (DEL, DUP) if the input SV file is a BED" }
  tx:
    type:
      - 'null'
      - type: enum
        name: tx
        symbols: ["RefSeq","ENSEMBL"]
    inputBinding:
      prefix: "-tx"
      position: 1
    doc: "Origin of the transcripts (RefSeq or ENSEMBL)"
  tx_file: { type: 'File?', inputBinding: { prefix: "-txFile", position: 1 }, doc: "file containing a list of preferred genes transcripts to be used in priority during the annotation (Preferred genes transcripts names should be tab or space separated)", "sbg:fileTypes": "BED, CSV, TSV, TXT" }

outputs:
  annotated_calls:
    type: 'File?'
    outputBinding:
      glob: '*.annotated.tsv'
      outputEval: $(inheritMetadata(self, inputs.sv_input_file))
  unannotated_calls:
    type: 'File?'
    outputBinding:
      glob: '*.unannotated.tsv'
      outputEval: $(inheritMetadata(self, inputs.sv_input_file))

$namespaces:
  sbg: https://sevenbridges.com
