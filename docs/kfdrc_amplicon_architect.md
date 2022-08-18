# Kids First DRC Amplicon Architect(AA) Workflow (WIP)
This is an extra-chromosomal DNA (ecDNA) pipeline designed to detect ecDNA candidates.
Briefly, chromosomal fragments in tumor cells break off to form circular DNA.
This DNA acts sort of like a plasmid/mtDNA in that it is not beholden to the usual rules of cell division and can have massive copy numbers in each cell.
The tools for this pipeline were derived from https://github.com/jluebeck/AmpliconArchitect, related publication: https://www.nature.com/articles/s41467-018-08200-y

## Workflow: `workflow/kfdrc_production_amplicon_architect.cwl`
This workflow allows you to run amplicon architect with various levels of "input readiness," geared toward working with outputs from the Kids First Somatic Workflow. It requires copy number info and the preferred method (by Kids First DRC) as input is from CNVkit. Given that the following perquisites must be fulfilled:
 - You have `.cns` output (_not_ `.call.cns`) readily available from a previous run
 - **OR** You have `.cnn` output from a previous run, allowing CNVkit to be run while skipping normal reference and annotation, saving time and money
 - You have the normal alignment (bam or cram) file. In this case, CNVkit will start from scratch
A rough run time estimate put this at about 4 hours and $1.40 in cost using spot instances on CAVATICA when CNVkit has already been run, to about 8 hours and $2.88 in cost. Run time varies if prepare step has candidates (and number of candidates) and availability of CNVkit inputs.
### Inputs (critical):
 - `aa_data_repo`: Reference tar ball obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/
 - `aa_data_ref_version`: enum, "GRCh38", "hg19", "GRCh37", "mm10", "GRCm38". Genome version in data repo to use. Should match `data_repo` downloaded.
`tumor_align_file`: "Tumor read alignment file. Can be cram or bam", secondaryFiles: [ '^.bai?', '.crai?' ] }
output_basename: { type: string, doc: "File name prefix for steps that require it" }
mosek_license_file: { type: File, doc: "This tool uses some software that requires a license file. You can get a personal or institutional one from https://www.mosek.com/license/request/." }
### Input (if cram input or .cnn not available):
 - `reference`: Reference fasta file
### Input (if `.cns` file available):
 - `cnvkit_cns`: .cns file from previous CNVkit run, if available. DO NOT USE .call.cns
### Input (if `.cnn` file available **and** `.cns` file _not_ available)
 - `cnvkit_cnn`:If running using an existing .cnn, supply here
### Input (if no CNVkit inputs available):
 - `normal_align_file`: Normal read alignment file. Can be cram or bam
 - `annotation_file`: refFlat.txt file, needed if cnv kit cnn not already built
 - `male_input_flag`: Boolean flag indicating if input is male or not
### Outputs:
 - `aa_cnv_seeds`: Bed file with candidate regions to search
 - `aa_summary`: Summary for all amplicons detected by AA
 - `aa_cycles`: Text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
 - `aa_graph`: A text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
 - `aa_sv_png`: PNG image file displaying the SV view of AA
 - `aa_classification_profiles`: Abstract classification of the amplicon
 - `aa_gene_list`: Genes present on amplicons with each classification

## Tools
### Dockerfiles:
 - `images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3`: Optional, to create a "raw" (not .call.cns) .cns file.
 - `jluebeck/prepareaa:v0.1203.10`: Contains all software required to run pre-processing and processing steps

### `tools/cnvkit_batch_only.cwl`
Modified version of somatic workflow `tools/cnvkit_batch.cwl` to truly run the batch step in a limited capacity.
Other CN callers can be used, but since KF uses this natively, but does not output `.cns` calls, it is necessary to run on legacy samples until this workflow is part of production.
#### Inputs (critical):
 - `input_sample`: tumor bam file
 - `run_mode`: enum with choices: ["wgs", "hybrid", "amplicon"], _currently AA works only with WGS_
 - `male_input_flag`: boolean, set to true if reference is male
 - `output_basename`: string
 If `.cnn` file available:
 - `cnvkit_cnn`: `.cnn` file from a previous run 
 If `.cnn` is unavailable:
 - `input_control`: input normal file
 - `reference`: fasta file and index
 - `annotation_file`: refFlat.txt file
#### Outputs (critical):
 - `output_cns`: `.cns` file needed for AA

### `tools/cns_to_aa_bed.cwl`
Converts `.cns` file to bed, adding "crude ploidy estimates".
Again, it is recommended that `.cns` be used, not `.call.cns`
#### Inputs:
 - `input_cns`: Raw `.cns` output
#### Outputs:
 - `aa_cns_bed`: Outputs two bed files. May restrict to one later once I understand better
   - `_ESTIMATED_PLOIDY_CORRECTED_CN.bed`: currently used for subsequent steps
   - `_uncorr_CN.bed`: currently not used, keep until I understand better

**Note, while the tools below can use cram input, processing is slower by about 16%**
### `tools/prepare_aa.cwl`
This tool uses the bed file regions to analyze input bam/cram file for amplicon candidates.
It is also capable of performing the functionality of the next two steps if the flags are given
#### Inputs (critical):
 - `data_repo`: Reference tar ball obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/
 - `data_ref_version`: enum, "GRCh38", "hg19", "GRCh37", "mm10", "GRCm38". Genome version in data repo to use. Should match `data_repo` downloaded.
 - `sorted_bam`: Input bam/cram
 - `sample`: Output prefix name
 - `cnv_bed`: Converted CNVkit cns-to-bed file
 - `ref_cache`: _For cram input_, provide tar ball of output of running misc/seq_cache_populate.pl from samtools on the reference fasta
#### Additional inputs (`run_aa` flag given):
 - `run_aa` Flag to run amplicon architect right after prepare_aa - set to True to do so
 - `mosek_license_file`: This tool uses some software that requires a license file. You can get a personal or institutional one from https://www.mosek.com/license/request/
#### Additional inputs (`run_ac` flag given):
 - `run_aa` Flag to run amplicon classifier right after amplicon architect - set to True to do so
#### Outputs:
 - `aa_cnv_seeds`: Bed files with seed intervals to search for viable ecDNA amplicons
 - `coverage_stats`: Coverage stats file to provide for running AA
#### Additional outputs (`run_aa` flag given):
 - `summary`: summary for all amplicons detected by AA
 - `cycles`: text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
 - `graph`: A text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
 - `sv_pdf`: PDF image file displaying the SV view of AA
#### Additional outputs (`run_ac` flag given):
 - `amplicon_classification_profiles`: abstract classification of the amplicon
 - `gene_list`: genes present on amplicons with each classification
 - `ecDNA_cts`: not yet defined

### `tools/amplicon_architect.cwl`
This ties up all the prep work and outputs ecDNA candidates based on the seed regions provided
#### Inputs (critical):
 - `data_repo`: Reference tar ball obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/
 - `data_ref_version`: enum, "GRCh38", "hg19", "GRCh37", "mm10", "GRCm38". Genome version in data repo to use. Should match `data_repo` downloaded.
 - `bam`: Input bam/cram
 - `downsample`: Recommended for anaylsis, recommended value be 10
 - `output_basename`: output file name prefix
 - `aa_bed`: Bed file from prepare-aa
 - `coverage_stats`: Coverage stats file generated by prepare-aa
 - `ref_cache`: _For cram input_, provide tar ball of output of running misc/seq_cache_populate.pl from samtools on the reference fasta
 - `mosek_license_file`: This tool uses some software that requires a license file. You can get a personal or institutional one from https://www.mosek.com/license/request/
 ### Outputs:
 - `summary`: summary for all amplicons detected by AA
 - `cycles`: text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
 - `graph`: A text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
 - `sv_pdf`: PDF image file displaying the SV view of AA
 - `sv_png`: PNG image file displaying the SV view of AA

### `tools/aa_classifier.cwl`
This tool classifies Amplicon Architect output to detect types of focal amplifications present
#### Inputs (critical):
 - `data_repo`: Reference tar ball obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/
 - `data_ref_version`: enum, "GRCh38", "hg19", "GRCh37", "mm10", "GRCm38". Genome version in data repo to use. Should match `data_repo` downloaded.
 - `cycles`: text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
 - `graph`: A text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
#### Outputs:
 - `amplicon_classification_profiles`: abstract classification of the amplicon
 - `gene_list`: genes present on amplicons with each classification
 - `ecDNA_cts`: not yet defined
