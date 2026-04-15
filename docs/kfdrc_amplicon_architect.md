# Kids First DRC Amplicon Architect(AA) Workflow
This is an extra-chromosomal DNA (ecDNA) pipeline designed to detect ecDNA candidates.
Briefly, chromosomal fragments in tumor cells break off to form circular DNA.
This DNA acts sort of like a plasmid/mtDNA in that it is not beholden to the usual rules of cell division and can have massive copy numbers in each cell.
The tools for this pipeline were derived from https://github.com/AmpliconSuite/AmpliconSuite-pipeline.

## [AmpliconSuite-pipeline](../tools/ampliconsuite_pipeline.cwl)
The primary work element of the workflow is a CWL wrapper for [AmpliconSuite-pipeline.py](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/v1.5.2/AmpliconSuite-pipeline.py).
The script is itself an intelligent pipeline that enables all steps (alignment, CNV calling, seed interval detection) prior to running AmpliconArchitect, and invokes AmpliconClassifier.
This tool can reasonably be run as a standalone element to achieve all the same results that one would get from the workflow.

## [Workflow](../workflow/kfdrc_production_amplicon_architect.cwl)
As mentioned above the AmpliconSuite-pipeline can simply be run to get results, so why a workflow? Mostly it comes down to minor timecost considerations. They are as follows:
- While modern CNVkit can handle CRAMs, it takes considerably longer to process them than BAMs. Even when accounting for the time to convert the CRAMs to BAMs before CNVkit, processing CRAMs takes longer and costs 20-30% more.
- While AmpliconSuite-pipeline can run both CNVkit and AmpliconArchitect in the same run, the pipeline breaks these up into separate steps. The reason for this is twofold:
  - CNVkit is multithreaded, AmpliconArchitect is single threaded. To limit the number of idle CPUs, we switch to a smaller instance for AmpliconArchitect
  - Both of these steps can be quite long and are vulnerable to spot instance kills. Breaking these up into separate steps reduces that vulnerability and better leverages memoization
### Inputs (critical):
 - `aa_data_repo`: Reference tar ball obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/
 - `tumor_align_file`: "Tumor read alignment file. Can be BAM or CRAM and must include index (BAI/CRAI)"
 - `output_basename`: "File name prefix for steps that require it"
 - `mosek_license_file`: "This workflow uses some software that requires a license file. You can get a personal or institutional one from https://www.mosek.com/license/request/."
### Input (if CRAM input):
 - `reference`: Reference fasta file
### Input (if `.cns` file available):
 - `cnvkit_cns`: .cns file from previous CNVkit run, if available. DO NOT USE .call.cns
### Input (if no CNVkit inputs available):
 - `normal_align_file`: Normal read alignment file. Can be CRAM or BAM
### Outputs:
 - `aa_cnv_seeds`: Bed file with candidate regions to search
 - `aa_summary`: Summary for all amplicons detected by AA
 - `aa_cycles`: Text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
 - `aa_graph`: A text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
 - `aa_sv_png`: PNG image file displaying the SV view of AA
 - `aa_classification_profiles`: Abstract classification of the amplicon
 - `aa_gene_list`: Genes present on amplicons with each classification
