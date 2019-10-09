class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://sevenbridges.com'
id: brownm28/mb-controlfreec-troubleshoot/control-freec-11-6-sbg/0
baseCommand: []
inputs:
  - 'sbg:category': General
    id: GC_content_profile
    type: File?
    label: GC content profile
    doc: GC-content profile for a given window-size.
    'sbg:fileTypes': CNP
  - 'sbg:category': General
    'sbg:toolDefaultValue': 'FALSE'
    id: bed_graph_output
    type: boolean?
    label: Bed Graph output
    doc: >-
      Set "BedGraphOutput=TRUE" if you want an additional output in BedGraph
      format for the UCSC genome browser.
  - 'sbg:category': General
    'sbg:toolDefaultValue': '0.8'
    id: break_point_threshold
    type: float?
    label: Break point threshold
    doc: >-
      Positive value of threshold for segmentation of normalized profiles. The
      closer it is to zero, the more breakpoints will be called. Its recommended
      value is between 0.1 and 1.2.
  - 'sbg:category': General
    'sbg:toolDefaultValue': '2'
    id: break_point_type
    type:
      - 'null'
      - type: enum
        symbols:
          - '0'
          - '1'
          - '2'
          - '3'
          - '4'
        name: break_point_type
    label: Break point type
    doc: >-
      Desired behavior in the ambiguous regions (poly-N or low mappability
      regions between two different copy number values). 0: the "unknown" region
      is attached to the "known" region on the right  1: make a separate
      fragment of this “unknown” region and then attaches it to the left or to
      the right region choosing the longer one  2: make a separate fragment of
      this “unknown” region and then attaches it to the left or to the right
      region but the “ploidy” copy number has a priority  3: make a separate
      fragment of this “unknown” region and then attaches it to the left or to
      the right region choosing the longer one but this “known” region should
      make at least half-size of the “unknown” region  4: make a separate
      fragment of this “unknown” region and do not assign any copy number to
      this region at all
  - 'sbg:category': Target
    id: capture_regions
    type: File?
    label: Capture regions
    doc: >-
      Capture regions in .bed format; sorted .bed file should contain the
      following colomns: chr   0-based start   1-based end.
    'sbg:fileTypes': BED
  - 'sbg:category': File Input
    id: chr_len
    type: File
    label: Chromosomes length file
    doc: Chromosome length file in a tab-delimited format.
    'sbg:fileTypes': 'TXT, LEN, SIZES'
  - 'sbg:toolDefaultValue': '0.05'
    id: coeff_var
    type: float?
    label: Coefficient of variation
    doc: Coefficient of variation to evaluate necessary window size.
  - 'sbg:category': General
    'sbg:toolDefaultValue': '0'
    id: contamination
    type: float?
    label: Contamination
    doc: >-
      A priori known value of tumor sample contamiantion by normal cells. Set
      "contaminationAdjustment=TRUE" to correct for the contamination.
  - 'sbg:category': General
    'sbg:toolDefaultValue': 'FALSE'
    id: contamination_adjustment
    type: boolean?
    label: Contamination adjustment
    doc: >-
      Set TRUE to correct for contamination by normal cells. If "contamination"
      is not provided, it will automatically evaluate the level of
      contamination.
  - 'sbg:category': General
    'sbg:toolDefaultValue': >-
      3&4 (GC-content based normalization, WGS) or 1 (control-read-count-based
      normalization, WES)
    id: degree
    type: float?
    label: Degree
    doc: Degree of polynomial.
  - 'sbg:category': General
    'sbg:toolDefaultValue': 'WGS: 0, WES: 1'
    id: force_GC_content_normalization
    type:
      - 'null'
      - type: enum
        symbols:
          - '0'
          - '1'
          - '2'
        name: force_GC_content_normalization
    label: Force GC content normalization
    doc: >-
      Set 1 or 2 to correct the Read Count (RC) for GC-content bias and low
      mappability even when you have a control sample. 0: simply model "sample
      RC ~ Control RC"  1: normalize the sample and the control RC using
      GC-content and then calculate the ratio "Sample RC/contol RC" 2: model
      "sample RC ~ Control RC" bias, and then normalize for GC-content.
  - 'sbg:category': File Input
    id: gem_mappability_file
    type: File?
    label: GEM mappability file
    doc: Mappability file in GEM format.
    'sbg:fileTypes': GEM
  - 'sbg:category': Execution
    'sbg:toolDefaultValue': '1 - with GC-content, 0 - with a control dataset'
    id: intercept
    type: float?
    label: Intercept
    doc: Intercept of polynomial.
  - 'sbg:category': Control
    id: mate_copynumber_file_control
    type: File?
    label: Mate copy number file control
    doc: >-
      Raw copy number profile for a given window-size (higher priority than
      mateFile)  (don't need to provide a mateFile if mateCopyNumberFile is
      provided).
    'sbg:fileTypes': CPN
  - 'sbg:category': Sample
    id: mate_copynumber_file_sample
    type: File?
    label: Mate copy number file sample
    doc: >-
      Raw copy number profile for a given window-size (higher priority than
      mateFile)  (don't need to provide a mateFile if mateCopyNumberFile is
      provided).
    'sbg:fileTypes': CPN
  - 'sbg:category': Control
    id: mate_file_control
    type: File?
    label: Mate file Control
    doc: 'Mapped reads (can be single end reads, mate-pairs or paired-end reads).'
    'sbg:fileTypes': 'SAM, BAM, PILEUP, PILEUP.GZ'
  - 'sbg:category': Sample
    id: mate_file_sample
    type: File?
    label: Mate file Sample
    doc: 'Mapped reads (can be single end reads, mate-pairs or paired-end reads).'
    'sbg:fileTypes': 'SAM, BAM, PILEUP, PILEUP.GZ'
  - 'sbg:category': Control
    id: mate_orientation_control
    type:
      - 'null'
      - type: enum
        symbols:
          - '0'
          - RF
          - FR
          - FF
        name: mate_orientation_control
    label: Mate orientation control
    doc: >-
      Format of reads (in mateFile). 0 (for single ends), RF (Illumina
      mate-pairs),  FR (Illumina paired-ends), FF (SOLiD mate-pairs).
  - 'sbg:category': Sample
    id: mate_orientation_sample
    type:
      - 'null'
      - type: enum
        symbols:
          - '0'
          - RF
          - FR
          - FF
        name: mate_orientation_sample
    label: Mate orientation sample
    doc: >-
      Format of reads (in mateFile). 0 (for single ends), RF (Illumina
      mate-pairs),  FR (Illumina paired-ends), FF (SOLiD mate-pairs).
  - 'sbg:category': General
    'sbg:toolDefaultValue': 0.55 (change only if you run Control-FREEC on a bacterial genome)
    id: max_expected_GC
    type: float?
    label: Maximum expected GC
    doc: >-
      Maximal exptected value of the GC-content for the prior evaluation of
      "Read Count ~ GC-content" dependancy.
  - 'sbg:category': General
    'sbg:toolDefaultValue': '8'
    id: max_threads
    type: int?
    label: Maximum threads
    doc: Number of threads (multi-threading mode).
  - 'sbg:category': General
    'sbg:toolDefaultValue': 'WES: 3, WGS: 1'
    id: min_CNA_length
    type: int?
    label: Minimum CNA length
    doc: Minimal number of consecutive windows to call a CNA.
  - 'sbg:category': General
    'sbg:toolDefaultValue': 0.35 (change only if you run Control-FREEC on a bacterial genome)
    id: min_expected_GC
    type: float?
    label: Minimal exptected GC
    doc: >-
      Minimal expected value of the GC-content for the prior evaluation of "Read
      Count ~ GC-content" dependency.
  - 'sbg:category': General
    'sbg:toolDefaultValue': '0.85'
    id: min_map_per_w
    type: float?
    label: Minimum mappability per window
    doc: >-
      Only windows with fraction of mappable positions higher than or equal to
      this threshold will be considered  (if "gemMappabilityFile" is not
      provided, one uses the percentage of non-N letters per window).
  - 'sbg:category': General
    'sbg:toolDefaultValue': >-
      100 (meaning "do not look for subclones"). Suggested: 20 (or 0.2) for WGS
      and 30 (or 0.3) for WES.
    id: min_subclone_presence
    type: int?
    label: Minimal subclone presence
    doc: Detects subclones present in x% of cell population
  - 'sbg:category': Control
    id: mini_pileup_control
    type: File?
    label: Mini pileup Control
    doc: >-
      Mini pileup file created from the corresponding BAM file dring a previous
      run of Control-FREEC - providing this file will significantly speed up the
      whole process.
    'sbg:fileTypes': PILEUP
  - 'sbg:category': Sample
    id: mini_pileup_sample
    type: File?
    label: Mini pileup Sample
    doc: >-
      Mini pileup file created from the corresponding BAM file dring a previous
      run of Control-FREEC - providing this file will significantly speed up the
      whole process.
    'sbg:fileTypes': PILEUP
  - 'sbg:category': BAF
    'sbg:toolDefaultValue': '0'
    id: minimal_coverage_per_position
    type: int?
    label: Minimal coverage per position
    doc: Minimal read coverage for a position to be considered in BAF analysis.
  - 'sbg:category': BAF
    'sbg:toolDefaultValue': '0'
    id: minimal_quality_per_position
    type: int?
    label: Minimal quality per position
    doc: >-
      Minimal sequencing quality for a position to be considered in BAF
      analysis. Default: 0; using this option can slow down reading of pileup
      files.
  - 'sbg:category': General
    'sbg:toolDefaultValue': 'FALSE'
    id: noisy_data
    type: boolean?
    label: Noisy data
    doc: >-
      Set TRUE for target resequencing data (e.g., exome-seq) to avoid false
      positive predictions due to nonuniform capture.
  - 'sbg:category': General
    id: ploidy
    type: 'int[]'
    label: Ploidy
    doc: >-
      Genome ploidy. In case of doubt, you can set different values and
      Control-FREEC will select the one that explains most observed CNAs (eg.
      2,3,4)
  - 'sbg:category': General
    'sbg:toolDefaultValue': 'TRUE'
    id: print_NA
    type: boolean?
    label: Print NA
    doc: >-
      Set FALSE to avoid printing "-1" to the _ratio.txt files Useful for
      exome-seq or targeted sequencing data.
  - 'sbg:category': General
    'sbg:toolDefaultValue': 10; recommended value >=50 for for exome data
    id: read_cnt_threshold
    type: int?
    label: Read count threshold
    doc: >-
      Threshold on the minimal number of reads per window in the control sample
      Useful for exome-seq or targeted sequencing data.
  - 'sbg:category': BAF
    id: reference
    type: File
    label: Reference file
    doc: >-
      Reference file that will be divided as needed to separate chromosomes and
      contigs.
    'sbg:fileTypes': 'FA, FASTA'
  - 'sbg:category': General
    id: sex
    type:
      - 'null'
      - type: enum
        symbols:
          - XX
          - XY
        name: sex
    label: Sex
    doc: Sample sex.
  - 'sbg:category': BAF
    'sbg:toolDefaultValue': '0'
    id: shift_in_quality
    type: int?
    label: Shift in quality
    doc: 'Basis for Phred quality. Default: 0; usually 33 or 64; see fastq quality.'
  - 'sbg:category': BAF
    id: snp_file
    type: File?
    label: Known SNPs
    doc: Known SNPs.
    'sbg:fileTypes': 'TXT, VCF'
  - 'sbg:category': General
    id: step
    type: int?
    label: Step
    doc: 'Step (used only when "window" is specified - Ex: 10000).'
  - 'sbg:category': General
    'sbg:toolDefaultValue': '50000'
    id: telocentromeric
    type: int?
    label: Telocentromeric
    doc: >-
      Length of pre-telomeric and pre-centromeric regions: Control-FREEC will
      not output small CNAs and LOH found within these regions (they are likely
      to be false because of mappability/genome assembly issues) 50000 is OK for
      human/mouse genomes. Use smaller values for yeasts and flies.
  - 'sbg:category': Execution
    id: total_memory
    type: int?
    label: 'Total memory [MB]'
    doc: Total amount of memory in MB reserved on the instance.
  - 'sbg:category': General
    'sbg:toolDefaultValue': 'FALSE'
    id: unique_match
    type: boolean?
    label: Unique match
    doc: >-
      Use a mappability profile to correct read counts (in this case a
      mappability file must be provided with "gemMappabilityFile").
  - id: window
    type: int?
    label: Window
    doc: explicit window size (higher priority than coefficientOfVariation)
outputs:
  - id: GC_profile
    doc: GC profile output file.
    label: GC profile
    type: File?
    outputBinding:
      glob: |-
        ${
            return "GC_profile.targetedRegions.cnp"
        }
      outputEval: |
        ${
            return inheritMetadata(self, inputs.mate_file_sample)

        }
    'sbg:fileTypes': CNP
  - id: cnvs
    doc: File with coordinates of predicted copy number alterations.
    label: CNVs output
    type: File?
    outputBinding:
      glob: '*_CNVs'
      outputEval: |+
        ${
            return inheritMetadata(self, inputs.mate_file_sample)

        }

    'sbg:fileTypes': TXT
  - id: cnvs_pvalue
    doc: File with coordinates of predicted copy number alterations with p-values.
    label: CNVs with p-value
    type: File?
    outputBinding:
      glob: '*.p.value.txt'
      outputEval: |
        ${
            return inheritMetadata(self, inputs.mate_file_sample)

        }
    'sbg:fileTypes': TXT
  - id: config_script
    doc: Configuration script used for running.
    label: Configuration script used for running
    type: File?
    outputBinding:
      glob: config.txt
      outputEval: |
        ${
            return inheritMetadata(self, inputs.mate_file_sample)

        }
    'sbg:fileTypes': TXT
  - id: control_cpn
    doc: Control CPN file to be used in the future for more efficient computation.
    label: Control CPN
    type: File?
    outputBinding:
      glob: '*_control.cpn'
      outputEval: |-
        ${  if (inputs.mate_file_control){
            return inheritMetadata(self, inputs.mate_file_control)
        }
        }
    'sbg:fileTypes': CPN
  - id: control_pileup
    doc: >-
      Mini Pileup created for BAF calculation. It is used to speed up consequent
      runs with the same samples.
    label: Control Pileup
    type: File?
    outputBinding:
      glob: |-
        ${
            if (inputs.mate_file_control) {

                return inputs.mate_file_control.path.split('/').pop() + '_minipileup.pileup'

            }
        }
      outputEval: |-
        ${
            if (inputs.mate_file_control){
            return inheritMetadata(self, inputs.mate_file_control)
        }
        }
    'sbg:fileTypes': PILEUP
  - id: info_txt
    doc: Parsable file with information about FREEC run.
    label: Info TXT
    type: File?
    outputBinding:
      glob: '*_info.txt'
      outputEval: |
        ${
            return inheritMetadata(self, inputs.mate_file_sample)

        }
    'sbg:fileTypes': TXT
  - id: pngs
    doc: >-
      Visalized normalized copy number profile with predicted CNAs as well as
      BAF profile (if dbSNP is provided)
    label: Copy number profile
    type: 'File[]?'
    outputBinding:
      glob: '*png'
    'sbg:fileTypes': PNG
  - id: ratio
    doc: File with ratios and predicted copy number alterations for each window.
    label: Ratio
    type: File?
    outputBinding:
      glob: '*_ratio.txt'
      outputEval: |
        ${
            return inheritMetadata(self, inputs.mate_file_sample)

        }
    'sbg:fileTypes': TXT
  - id: ratio_BedGraph
    doc: >-
      File with ratios in BedGraph format for visualization in the UCSC genome
      browser
    label: Ratio BedGraph
    type: File?
    outputBinding:
      glob: '*.BedGraph'
      outputEval: |
        ${
            return inheritMetadata(self, inputs.mate_file_sample)

        }
    'sbg:fileTypes': BEDGRAPH
  - id: sample_BAF
    doc: File with B-allele frequencies for each possibly heterozygous SNP position
    label: BAF sample file
    type: File?
    outputBinding:
      glob: '*_BAF.txt'
      outputEval: |
        ${
            return inheritMetadata(self, inputs.mate_file_sample)

        }
    'sbg:fileTypes': TXT
  - id: sample_cpn
    doc: Sample CPN file to be used in the future for more efficient computation.
    label: Sample CPN
    type: File?
    outputBinding:
      glob: '*_sample.cpn'
      outputEval: |-
        ${
            return inheritMetadata(self, inputs.mate_file_sample)

        }
    'sbg:fileTypes': CPN
  - id: sample_pileup
    doc: >-
      Mini Pileup created for BAF calculation. It is used to speed up consequent
      runs with the same samples.
    label: Sample Pileup
    type: File?
    outputBinding:
      glob: |-
        ${
            if (inputs.mate_file_sample) {

                return inputs.mate_file_sample.path.split('/').pop() + '_minipileup.pileup'
            }
        }
      outputEval: |-
        ${
            return inheritMetadata(self, inputs.mate_file_sample)

        }
    'sbg:fileTypes': PILEUP
doc: >-
  Control-FREEC analyzes copy-number variants and allelic imbalances in exome
  and whole-genome DNA sequencing.


  This tool automatically computes, normalizes and segments copy number and beta
  allele frequency (BAF) profiles, then calls copy number alterations and LOH.
  [1]


  *A list of **all inputs and parameters** with corresponding descriptions can
  be found at the bottom of the page.*

  ### Common Use Cases


  * The **chrLenFile** input is required and can be found in Public Reference
  Files as **Homo\_sapiens\_assembly38.fasta.sizes**, **ucsc.hg19.fasta.sizes**
  and **human\_g1k\_v37\_decoy.fasta.sizes**.


  * The **ploidy** parameter is required. In case of doubt, different values can
  be set and Control-FREEC will select the one that explains the most observed
  CNVs.


  * Normal and control sample can be provided through two possible inputs:
       * **mateFile**, a file with mapped reads
       * **mateCopyNumberFile**, a raw copy number file created for both normal and control sample, provided through **mateFile** in a first run, and can be reused in the future runs for more efficient computation.



  * **A control (matched normal) sample is optional for whole genome sequencing
  data but mandatory for whole exome or targeted sequencing data.**


  * Similar to **mateCopyNumberFile**, a **Mini pileup Sample** and **Mini
  pileup Control** files can be created in the first run, if the **Known SNPs**
  file is provided. Consequently, by providing these files as inputs in future
  tasks, execution time will decrease significantly.


  * If a **mateFile** is specified, the **mateOrientation** parameter must be
  set.


  * In order to create a **BAF profile**, one of the following options must be
  implemented:
      * **mateFile** + **Known SNPs** 
      * **mateCopyNumberFile** + **mateFile** + **KnownSNPs**
      * **mateCopyNumberFile** + **miniPileup** + **KnownSNPs**

  ### Changes Introduced by Seven Bridges


  * Based on the input parameters, a config file is created in order to properly
  run Control-FREEC.


  ### Common Issues and Important Notes


  * **A control (matched normal) sample is optional for whole genome sequencing
  data but mandatory for whole exome or targeted sequencing data.**


  * A **gemMappabilityFile** can be used only in the mode without a control
  sample.


  * If a **mateFile** is specified, the **mateOrientation** parameter must be
  set.


  * Currently, there is an issue with creating a **BAF sample file** with the
  b37 notation. The genotypes for CNV regions are, however, created.



  ### Performance Benchmarking


  The instance set for this tool is the AWS c4.2xlarge instance with 8 vCPUs, 15
  GiB of RAM and 1 TB of EBS (disk space).

  |     BAM size in GB    | Type |  Instance  | Duration | Cost ($) |

  |:--------------------:|:----:|:----------:|:--------:|:--------:|

  |         2x12 (Normal-Tumor)         |  WES | c4.2xlarge |  1h 52m  |   
  0.8   |

  |   100 (Tumor-only)   |  WGS | c4.2xlarge |  17h 43m |     7    |

  | 2x100 (Normal-Tumor) |  WGS | c4.2xlarge |   1d 8h  |    13    |

  |   100 (Tumor-only)   |  WGS | c4.8xlarge |  6h 30m |     10    |

  | 2x100 (Normal-Tumor) |  WGS | c4.8xlarge |   11h  |    18    |


  An instance with more resources can be obtained by providing inputs for
  **Maximum threads** and **Total memory [MB]**.


  *Cost can be significantly reduced by using **spot instances**. Visit the
  [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances)
  for more details.*  


  ###References

  [1] [Control-FREEC: Prediction of copy number alterations and loss of
  heterozygosity using deep-sequencing
  data](http://boevalab.com/FREEC/tutorial.html#install)
label: Control-FREEC 11.6
arguments:
  - position: 0
    shellQuote: false
    valueFrom: |-
      ${
          // script for splitting the genome fasta into chromosomes
          return 'python split_fasta.py ' + inputs.reference.path

      }
  - position: 1
    shellQuote: false
    valueFrom: '&&'
  - position: 2
    prefix: ''
    shellQuote: false
    valueFrom: /opt/controlfreec/FREEC/src/freec
  - position: 3
    shellQuote: false
    valueFrom: '-conf'
  - position: 4
    shellQuote: false
    valueFrom: config.txt
  - position: 5
    shellQuote: false
    valueFrom: '&&'
  - position: 6
    shellQuote: false
    valueFrom: |-
      ${




          if (inputs.mate_file_sample) {
              filepath = inputs.mate_file_sample.path
              filename = filepath.split("/").pop()
          } else {
              filepath = inputs.mate_copynumber_file_sample.path
              filename = filepath.split("/").pop()
          }

          CNVs = filename + "_CNVs"
          ratio = filename + "_ratio" + ".txt"


          return "cat assess_significance.R | R --slave --args " + CNVs + " " + ratio
      }
  - position: 7
    shellQuote: false
    valueFrom: '&&'
  - position: 8
    shellQuote: false
    valueFrom: |-
      ${
          return "line=$(cat *info.txt | grep Output_Ploidy | sed -E 's/.+([0-9]+)/\\1/')"
      }
  - position: 9
    shellQuote: false
    valueFrom: '&&'
  - position: 10
    shellQuote: false
    valueFrom: |-
      ${
          return "cat makeGraph.R | R --slave --args"
      }
  - position: 11
    shellQuote: false
    valueFrom: $line
  - position: 12
    shellQuote: false
    valueFrom: |-
      ${
          if (inputs.mate_file_sample) {
              filepath = inputs.mate_file_sample.path
              filename = filepath.split("/").pop()
          } else {
              filepath = inputs.mate_copynumber_file_sample.path
              filename = filepath.split("/").pop()
          }

          ratio = filename + "_ratio" + ".txt"

          return ratio
      }
  - position: 13
    shellQuote: false
    valueFrom: |-
      ${

          if (inputs.snp_file) {

              sufix = "_BAF"
              sufix_ext = ".txt"

              if (inputs.mate_file_sample) {
                  filepath = inputs.mate_file_sample.path
                  filename = filepath.split("/").pop()
              } else {
                  filepath = inputs.mate_copynumber_file_sample.path
                  filename = filepath.split("/").pop()
              }


              new_filename = filename + sufix + sufix_ext

              return new_filename
          }
      }
  - position: 114
    shellQuote: false
    valueFrom: |-
      ${ //conversion of file names

          if (inputs.mate_file_control) {
              if (inputs.mate_file_control.path.split('.').pop() != 'pileup') {
                  com = ''
                  com += '&& mv sample.pileup '

              }
          }


      }
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: |-
      ${
          if (inputs.total_memory) {
              return inputs.total_memory
          } else {
              return 15000
          }
      }
    coresMin: |-
      ${
          if (inputs.max_threads) {
              return inputs.max_threads
          } else {
              return 8
          }
      }
  - class: DockerRequirement
    dockerImageId: caf6947244fa
    dockerPull: 'images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1'
  - class: InitialWorkDirRequirement
    listing:
      - entryname: config.txt
        entry: |-
          ${

              // The function returns the concatenated line for config file
              function makeline(content, p1, p2) {
                  if (p2 != null) {
                      if (p2.path != null) {
                          p2 = p2.path
                      }
                      content = content.concat(p1)
                      content = content.concat(" = ")
                      content = content.concat(p2)
                      content = content.concat("\n")
                  }
                  return content

              }

              // General section
              content = "[general]\n\n"
              content = makeline(content, "BedGraphOutput", inputs.bed_graph_output)
              content = content.concat("bedtools = /opt/bedtools2/bin/bedtools\n")
              content = makeline(content, "chrLenFile", inputs.chr_len)
              content = makeline(content, "breakPointThreshold", inputs.break_point_threshold)
              content = makeline(content, "breakPointType", inputs.break_point_type)
              content = makeline(content, "chrFiles", ".")
              content = makeline(content, "coefficientOfVariation", inputs.coeff_var)

              if (inputs.capture_regions) {
                  content = content.concat("window = 0\n")
              } else {
                  content = makeline(content, "window", inputs.window)
              }

              content = makeline(content, "contamination", inputs.contamination)
              content = makeline(content, "contaminationAdjustment", inputs.contamination_adjustment)
              content = makeline(content, "degree", inputs.degree)
              content = makeline(content, "forceGCcontentNormalization", inputs.force_GC_content_normalization)
              content = makeline(content, "GCcontentProfile", inputs.GC_content_profile)
              content = makeline(content, "gemMappabilityFile", inputs.gem_mappability_file)
              content = makeline(content, "intercept", inputs.intercept)
              content = makeline(content, "minCNAlength", inputs.min_CNA_length)
              content = makeline(content, "minMappabilityPerWindow", inputs.min_map_per_w)
              content = makeline(content, "minExpectedGC", inputs.min_expected_GC)
              content = makeline(content, "maxExpectedGC", inputs.max_expected_GC)
              content = makeline(content, "minimalSubclonePresence", inputs.min_subclone_presence)
              if (inputs.max_threads) {
                  content = makeline(content, "maxThreads", inputs.max_threads)
                  content = makeline(content, "SambambaThreads", inputs.max_threads)
              } else {
                  content = content.concat("maxThreads = 8\n")
                  content = content.concat("SambambaThreads = 8\n")
              }
              content = makeline(content, "noisyData", inputs.noisy_data)
              content = makeline(content, "ploidy", inputs.ploidy.toString())
              content = makeline(content, "printNA", inputs.print_NA)
              content = makeline(content, "readCountThreshold", inputs.read_cnt_threshold)
              content = content.concat("sambamba = /opt/sambamba_0.5.9/sambamba_v0.5.9\n")

              content = content.concat("samtools = /opt/samtools-1.3.1/samtools\n")
              content = makeline(content, "sex", inputs.sex)
              content = makeline(content, "step", inputs.step)
              content = makeline(content, "telocentromeric", inputs.telocentromeric)
              content = makeline(content, "uniqueMatch", inputs.unique_match)


              // Sample section

              content = content.concat("\n[sample]\n\n")
              content = makeline(content, "mateFile", inputs.mate_file_sample)
              content = makeline(content, "mateCopyNumberFile", inputs.mate_copynumber_file_sample)
              content = makeline(content, "miniPileup", inputs.mini_pileup_sample)
              if (inputs.mate_file_sample) {
                  if (inputs.mate_file_sample.path.split('.').pop() == "gz") {
                      content = makeline(content, "inputFormat", inputs.mate_file_sample.path.split('.').slice(-2, -1)[0])
                  } else {
                      content = makeline(content, "inputFormat", inputs.mate_file_sample.path.split('.').pop())
                  }
                  content = makeline(content, "mateOrientation", inputs.mate_orientation_sample)
              }


              // Control section

              content = content.concat("\n[control]\n\n")
              content = makeline(content, "mateFile", inputs.mate_file_control)
              content = makeline(content, "mateCopyNumberFile", inputs.mate_copynumber_file_control)
              content = makeline(content, "miniPileup", inputs.mini_pileup_control)
              if (inputs.mate_file_control) {
                  content = makeline(content, "inputFormat", inputs.mate_file_control.path.split('.').pop())
                  content = makeline(content, "mateOrientation", inputs.mate_orientation_sample)
              }




              // BAF section

              content = content.concat("\n[BAF]\n\n")
              content = makeline(content, "minimalCoveragePerPosition", inputs.minimal_coverage_per_position)
              content = makeline(content, "minimalQualityPerPosition", inputs.minimal_quality_per_position)
              content = makeline(content, "shiftInQuality", inputs.shift_in_quality)
              if (inputs.snp_file) {
                  content = makeline(content, "SNPfile", inputs.snp_file)
                  if (inputs.mate_file_sample) {
                      if ((inputs.mate_file_sample.path.split('.').pop().toUpperCase() != 'PILEUP') &&
                          (inputs.mate_file_sample.path.split('.').slice(-2, -1)[0].toUpperCase() != 'PILEUP')) {
                          content = makeline(content, "makePileup", inputs.snp_file)
                          content = makeline(content, "fastaFile", inputs.reference)
                      }
                  }
              }

              // Target section

              content = content.concat("\n[target]\n\n")
              content = makeline(content, "captureRegions", inputs.capture_regions)

              return content
          }
        writable: false
      - entryname: split_fasta.py
        entry: |-
          import sys

          with open(sys.argv[1], "r") as f:
              fasta = f.readlines()

          ref_lines = {}
          for i in range(0, len(fasta)):
              if fasta[i][0] == ">":
                  chrom = fasta[i].split()[0].split(">")[1]
                  print("Reading chromosome: " + chrom)
                  ref_lines[chrom] = [fasta[i]]
              else:
                  ref_lines[chrom].append(fasta[i])

          for chromosome, lines in ref_lines.items():
              print("Creating " + chromosome + ".fasta")
              with open(chromosome + ".fasta", "w") as chr_fasta:
                  for line in lines:
                      chr_fasta.write(line)
        writable: false
      - entryname: assess_significance.R
        entry: "#!/usr/bin/env Rscript\n\nlibrary(rtracklayer)\n\nargs <- commandArgs()\n\ndataTable <-read.table(args[5], header=TRUE);\nratio<-data.frame(dataTable)\n\ndataTable <- read.table(args[4], header=FALSE)\ncnvs<- data.frame(dataTable) \n\nratio$Ratio[which(ratio$Ratio==-1)]=NA\n\ncnvs.bed=GRanges(cnvs[,1],IRanges(cnvs[,2],cnvs[,3]))  \nratio.bed=GRanges(ratio$Chromosome,IRanges(ratio$Start,ratio$Start),score=ratio$Ratio)\n\noverlaps <- subsetByOverlaps(ratio.bed,cnvs.bed)\nnormals <- setdiff(ratio.bed,cnvs.bed)\nnormals <- subsetByOverlaps(ratio.bed,normals)\n\n#mu <- mean(score(normals),na.rm=TRUE)\n#sigma<- sd(score(normals),na.rm=TRUE)\n\n#hist(score(normals),n=500,xlim=c(0,2))\n#hist(log(score(normals)),n=500,xlim=c(-1,1))\n\n#shapiro.test(score(normals)[which(!is.na(score(normals)))][5001:10000])\n#qqnorm (score(normals)[which(!is.na(score(normals)))],ylim=(c(0,10)))\n#qqline(score(normals)[which(!is.na(score(normals)))], col = 2)\n\n#shapiro.test(log(score(normals))[which(!is.na(score(normals)))][5001:10000])\n#qqnorm (log(score(normals))[which(!is.na(score(normals)))],ylim=(c(-6,10)))\n#qqline(log(score(normals))[which(!is.na(score(normals)))], col = 2)\n\nnumberOfCol=length(cnvs)\n\nfor (i in c(1:length(cnvs[,1]))) {\n  values <- score(subsetByOverlaps(ratio.bed,cnvs.bed[i]))\n  #wilcox.test(values,mu=mu)\n  W <- function(values,normals){resultw <- try(wilcox.test(values,score(normals)), silent = TRUE)\n\tif(class(resultw)==\"try-error\") return(list(\"statistic\"=NA,\"parameter\"=NA,\"p.value\"=NA,\"null.value\"=NA,\"alternative\"=NA,\"method\"=NA,\"data.name\"=NA)) else resultw}\n  KS <- function(values,normals){resultks <- try(ks.test(values,score(normals)), silent = TRUE)\n\tif(class(resultks)==\"try-error\") return(list(\"statistic\"=NA,\"p.value\"=NA,\"alternative\"=NA,\"method\"=NA,\"data.name\"=NA)) else resultks}\n  #resultks <- try(KS <- ks.test(values,score(normals)), silent = TRUE)\n  #\tif(class(resultks)==\"try-error\") NA) else resultks\n  cnvs[i,numberOfCol+1]=W(values,normals)$p.value\n  cnvs[i,numberOfCol+2]=KS(values,normals)$p.value\n  }\n\nif (numberOfCol==5) {\n  names(cnvs)=c(\"chr\",\"start\",\"end\",\"copy number\",\"status\",\"WilcoxonRankSumTestPvalue\",\"KolmogorovSmirnovPvalue\")  \n}\nif (numberOfCol==7) {\n  names(cnvs)=c(\"chr\",\"start\",\"end\",\"copy number\",\"status\",\"genotype\",\"uncertainty\",\"WilcoxonRankSumTestPvalue\",\"KolmogorovSmirnovPvalue\")  \n}\nif (numberOfCol==9) {\n  names(cnvs)=c(\"chr\",\"start\",\"end\",\"copy number\",\"status\",\"genotype\",\"uncertainty\",\"somatic/germline\",\"precentageOfGermline\",\"WilcoxonRankSumTestPvalue\",\"KolmogorovSmirnovPvalue\")  \n}\nwrite.table(cnvs, file=paste(args[4],\".p.value.txt\",sep=\"\"),sep=\"\\t\",quote=F,row.names=F)"
        writable: false
      - entryname: makeGraph.R
        entry: "#!/usr/bin/env Rscript\n\nargs <- commandArgs()\n\ndataTable <-read.table(args[5], header=TRUE);\n\nratio<-data.frame(dataTable)\nploidy <- type.convert(args[4])\n\n\npng(filename = paste(args[5],\".log2.png\",sep = \"\"), width = 1180, height = 1180,\n    units = \"px\", pointsize = 20, bg = \"white\", res = NA)\nplot(1:10)\nop <- par(mfrow = c(5,5))\n\nfor (i in c(1:22,'X','Y')) {\n\ttt <- which(ratio$Chromosome==i)\n\tif (length(tt)>0) {\n\t plot(ratio$Start[tt],log2(ratio$Ratio[tt]),xlab = paste (\"position, chr\",i),ylab = \"normalized copy number profile (log2)\",pch = \".\",col = colors()[88])\n\t tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )\n\t points(ratio$Start[tt],log2(ratio$Ratio[tt]),pch = \".\",col = colors()[136])\n\t\n\t\n\ttt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)\n\t points(ratio$Start[tt],log2(ratio$Ratio[tt]),pch = \".\",col = colors()[461])\n\t tt <- which(ratio$Chromosome==i)\n\t \n\t #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:\n\t #points(ratio$Start[tt],log2(ratio$CopyNumber[tt]/ploidy), pch = \".\", col = colors()[24],cex=4)\n\t \n\t}\n\ttt <- which(ratio$Chromosome==i)\n\t\n\t#UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:\n\tpoints(ratio$Start[tt],log2(ratio$MedianRatio[tt]), pch = \".\", col = colors()[463],cex=4)\n\t\n}\n\ndev.off()\n\n\npng(filename = paste(args[5],\".png\",sep = \"\"), width = 1180, height = 1180,\n    units = \"px\", pointsize = 20, bg = \"white\", res = NA)\nplot(1:10)\nop <- par(mfrow = c(5,5))\n\nmaxLevelToPlot <- 3\nfor (i in c(1:length(ratio$Ratio))) {\n\tif (ratio$Ratio[i]>maxLevelToPlot) {\n\t\tratio$Ratio[i]=maxLevelToPlot;\n\t}\n}\n\n\nfor (i in c(1:22,'X','Y')) {\n\ttt <- which(ratio$Chromosome==i)\n\tif (length(tt)>0) {\n\t plot(ratio$Start[tt],ratio$Ratio[tt]*ploidy,ylim = c(0,maxLevelToPlot*ploidy),xlab = paste (\"position, chr\",i),ylab = \"normalized copy number profile\",pch = \".\",col = colors()[88])\n\t tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )\n\t points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = \".\",col = colors()[136])\n\t\n\ttt <- which(ratio$Chromosome==i  & ratio$Ratio==maxLevelToPlot & ratio$CopyNumber>ploidy)\t\n\tpoints(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = \".\",col = colors()[136],cex=4)\n\t \n\ttt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)\n\t points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = \".\",col = colors()[461])\n\t tt <- which(ratio$Chromosome==i)\n\t \n\t #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:\n\t #points(ratio$Start[tt],ratio$CopyNumber[tt], pch = \".\", col = colors()[24],cex=4)\n\t \n\t}\n\ttt <- which(ratio$Chromosome==i)\n\t\n\t#UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:\n\tpoints(ratio$Start[tt],ratio$MedianRatio[tt]*ploidy, pch = \".\", col = colors()[463],cex=4)\n\t\n}\n\ndev.off()\n\n\n\n\nif (length(args)>5) {\n\tdataTable <-read.table(args[6], header=TRUE);\n\tBAF<-data.frame(dataTable)\n\n\tpng(filename = paste(args[6],\".png\",sep = \"\"), width = 1180, height = 1180,\n\t    units = \"px\", pointsize = 20, bg = \"white\", res = NA)\n\tplot(1:10)\n\top <- par(mfrow = c(5,5))\n\n\tfor (i in c(1:22,'X','Y')) {\n\t    tt <- which(BAF$Chromosome==i)\n\t    if (length(tt)>0){\n\t\tlBAF <-BAF[tt,]\n\t\tplot(lBAF$Position,lBAF$BAF,ylim = c(-0.1,1.1),xlab = paste (\"position, chr\",i),ylab = \"BAF\",pch = \".\",col = colors()[1])\n\n\t\ttt <- which(lBAF$A==0.5)\t\t\n\t\tpoints(lBAF$Position[tt],lBAF$BAF[tt],pch = \".\",col = colors()[92])\n\t\ttt <- which(lBAF$A!=0.5 & lBAF$A>=0)\n\t\tpoints(lBAF$Position[tt],lBAF$BAF[tt],pch = \".\",col = colors()[62])\n\t\ttt <- 1\n\t\tpres <- 1\n\n\t\tif (length(lBAF$A)>4) {\n\t\t\tfor (j in c(2:(length(lBAF$A)-pres-1))) {\n\t\t\t\tif (lBAF$A[j]==lBAF$A[j+pres]) {\t\n\t\t\t\t\ttt[length(tt)+1] <- j \n\t\t\t\t}\n\t\t\t}\n\t\t\tpoints(lBAF$Position[tt],lBAF$A[tt],pch = \".\",col = colors()[24],cex=4)\n\t\t\tpoints(lBAF$Position[tt],lBAF$B[tt],pch = \".\",col = colors()[24],cex=4)\t\n\t\t}\n\n\t\ttt <- 1\n\t\tpres <- 1\n\t\tif (length(lBAF$FittedA)>4) {\n\t\t\tfor (j in c(2:(length(lBAF$FittedA)-pres-1))) {\n\t\t\t\tif (lBAF$FittedA[j]==lBAF$FittedA[j+pres]) {\t\n\t\t\t\t\ttt[length(tt)+1] <- j \n\t\t\t\t}\n\t\t\t}\n\t\t\tpoints(lBAF$Position[tt],lBAF$FittedA[tt],pch = \".\",col = colors()[463],cex=4)\n\t\t\tpoints(lBAF$Position[tt],lBAF$FittedB[tt],pch = \".\",col = colors()[463],cex=4)\t\n\t\t}\n\n\t   }\n\n\t}\n\tdev.off()\n\n}"
        writable: false
  - class: InlineJavascriptRequirement
    expressionLib:
      - |-
        var updateMetadata = function(file, key, value) {
            file['metadata'][key] = value;
            return file;
        };


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

        var toArray = function(file) {
            return [].concat(file);
        };

        var groupBy = function(files, key) {
            var groupedFiles = [];
            var tempDict = {};
            for (var i = 0; i < files.length; i++) {
                var value = files[i]['metadata'][key];
                if (value in tempDict)
                    tempDict[value].push(files[i]);
                else tempDict[value] = [files[i]];
            }
            for (var key in tempDict) {
                groupedFiles.push(tempDict[key]);
            }
            return groupedFiles;
        };

        var orderBy = function(files, key, order) {
            var compareFunction = function(a, b) {
                if (a['metadata'][key].constructor === Number) {
                    return a['metadata'][key] - b['metadata'][key];
                } else {
                    var nameA = a['metadata'][key].toUpperCase();
                    var nameB = b['metadata'][key].toUpperCase();
                    if (nameA < nameB) {
                        return -1;
                    }
                    if (nameA > nameB) {
                        return 1;
                    }
                    return 0;
                }
            };

            files = files.sort(compareFunction);
            if (order == undefined || order == "asc")
                return files;
            else
                return files.reverse();
        };
successCodes:
  - 0
temporaryFailCodes:
  - 1
'sbg:categories':
  - Copy-Number-Analysis
'sbg:cmdPreview': >-
  python split_fasta.py pat/sbgpro/path/to/reference.ext &&
  /opt/controlfreec/FREEC-11.5/src/freec -conf config.txt && cat
  assess_significance.R | R --slave --args chr_19.noDup0.pileup.gz_CNVs
  chr_19.noDup0.pileup.gz_ratio.txt && line=$(cat *info.txt | grep Output_Ploidy
  | sed -E 's/.+([0-9]+)/\1/') && cat makeGraph.R | R --slave --args $line
  chr_19.noDup0.pileup.gz_ratio.txt chr_19.noDup0.pileup.gz_BAF.txt
'sbg:image_url': null
'sbg:license': GNU General Public License v3.0 only
'sbg:links':
  - id: 'http://bioinfo-out.curie.fr/projects/freec/'
    label: Homepage
  - id: 'https://github.com/BoevaLab/FREEC'
    label: GitHub
  - id: >-
      http://bioinformatics.oxfordjournals.org/content/early/2011/12/05/bioinformatics.btr670
    label: Paper
'sbg:toolAuthor': Bioinformatics Laboratory of Institut Curie
'sbg:toolkit': Control-FREEC
'sbg:toolkitVersion': '11.6'
'sbg:projectName': MB ControlFREEC Troubleshoot
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedBy': brownm28
    'sbg:modifiedOn': 1569186036
    'sbg:revisionNotes': Copy of bogdang/controlfreec/control-freec-11-6/1
'sbg:appVersion':
  - v1.0
'sbg:id': brownm28/mb-controlfreec-troubleshoot/control-freec-11-6-sbg/0
'sbg:revision': 0
'sbg:revisionNotes': Copy of bogdang/controlfreec/control-freec-11-6/1
'sbg:modifiedOn': 1569186036
'sbg:modifiedBy': brownm28
'sbg:createdOn': 1569186036
'sbg:createdBy': brownm28
'sbg:project': brownm28/mb-controlfreec-troubleshoot
'sbg:sbgMaintained': false
'sbg:validationErrors': []
'sbg:contributors':
  - brownm28
'sbg:latestRevision': 0
'sbg:publisher': sbg
'sbg:content_hash': a8c609371987909ebbcad810349a09839fd48c2a7d4481c277e28cce897ed68cd
'sbg:copyOf': bogdang/controlfreec/control-freec-11-6/1
