cwlVersion: v1.2
class: Workflow
id: kfdrc-somatic-variant-workflow
label: Kids First DRC Somatic Variant Workflow
doc: |
  # Kids First DRC Somatic Variant Workflow

  <p align="center">
    <img src="https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png">
  </p>

  This repository contains the Kids First Data Resource Center (DRC) Somatic Variant Workflow, which includes somatic variant (SNV), copy number variation (CNV), and structural variant (SV) calls.
  The [somatic variant workflow](./workflow/kfdrc-somatic-variant-workflow.cwl) takes aligned cram input and performs somatic variant calling using Strelka2, Mutect2, Lancet, and VarDict Java, CNV estimation using Control-FREEC, CNVkit, and GATK, and SV calls using Manta.
  For whole genome sequencing data, the workflow will also predict extra chromosomal DNA (ecDNA) using AmpliconArchitect
  Somatic variant call results are annotated with hotspots, assigned population frequencies using gnomAD AF, calculated gene models using Variant Effect Predictor (VEP), then added an additional MAF output using a modified version of Memorial Sloan Kettering Cancer Center's (MSKCC) vcf2maf.
  See [annotation workflow doc](./docs/kfdrc_annotation_wf.md) for more details on annotation.

  If you would like to run this workflow using the CAVATICA public app, a basic primer on running public apps can be found [here](https://www.notion.so/d3b/Starting-From-Scratch-Running-Cavatica-af5ebb78c38a4f3190e32e67b4ce12bb).
  Alternatively, if you'd like to run it locally using `cwltool`, a basic primer on that can be found [here](https://www.notion.so/d3b/Starting-From-Scratch-Running-CWLtool-b8dbbde2dc7742e4aff290b0a878344d) and combined with app-specific info from the readme below.
  This workflow is the current production workflow, equivalent to this [CAVATICA public app](https://cavatica.sbgenomics.com/public/apps#cavatica/apps-publisher/kfdrc-somatic-variant-workflow).


  ## Somatic Variant Workflow Callers

  | Caller                                                  | [CNV](#cnv-estimators) | [SNV](#snv-callers) | [SV](#sv-callers) | [ecDNA](#ecdna-prediction) |
  |---------------------------------------------------------|:----------------------:|:-------------------:|:-----------------:|:--------------------------:|
  | [AmpliconArchitect](./docs/kfdrc_amplicon_architect.md) |                        |                     |                   |              x             |
  | CNVkit                                                  |            x           |                     |                   |                            |
  | Control-FREEC                                           |            x           |                     |                   |                            |
  | [GATK CNV](./docs/kfdrc_gatk_cnv_somatic_pair_wf.md)    |            x           |                     |                   |                            |
  | [Lancet](./docs/kfdrc_lancet_sub_wf.md)                 |                        |          x          |                   |                            |
  | Manta                                                   |                        |                     |         x         |                            |
  | [Mutect2](./docs/kfdrc_mutect2_sub_wf.md)               |                        |          x          |                   |                            |
  | [Strelka2](./docs/kfdrc_strelka2_subworkflow.md)        |                        |          x          |                   |                            |
  | THeTa2                                                  |                        |                     |         x         |                            |
  | [VarDict](./docs/kfdrc_vardict_sub_wf.md)               |                        |          x          |                   |                            |

  ### SNV Callers

  - [Strelka2](https://github.com/Illumina/strelka/tree/v2.9.3) `v2.9.3`, from Illumina, calls SNVs and insertions/deletions (INDEL)
    - See the [subworkflow doc](./docs/kfdrc_strelka2_subworkflow.md) for more information
  - [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360036730411-Mutect2) `v4.1.1.0`, from the Broad Institute, calls SNVs, multi-nucleotide variants (MNV: equal length substitutions over 1 base long), and INDELs
    - See the [subworkflow doc](./docs/kfdrc_mutect2_sub_wf.md) for more information
  - [Lancet](https://github.com/nygenome/lancet/releases/tag/v1.0.7) `v1.0.7`, from the New York Genome Center (NYGC), calls SNVs, MNVs, and INDELs
    - See the [subworkflow doc](./docs/kfdrc_lancet_sub_wf.md) for more information
  - [VarDict Java](https://github.com/AstraZeneca-NGS/VarDictJava/tree/1.7.0) `v1.7.0`, from AstraZeneca, calls SNVs, MNVs, INDELs and more
    - See the [subworkflow doc](./docs/kfdrc_vardict_sub_wf.md) for more information

  #### SNV Consensus Calling

  Each caller has a different approach to variant calling, and together one can glean confident results. To synthesize the results from these four callers we have developed a [consensus calling workflow](./workflow/kfdrc_consensus_calling.cwl). For details, please refer to the [workflow documentation](./docs/kfdrc-consensus-calling.md).

  #### SNV Annotation

  Somatic variant call results are annotated with hotspots, assigned population frequencies using gnomAD AF, calculated gene models using VEP, then added an additional MAF output using a modified version of MSKCC vcf2maf.

  The somatic annotation workflow is included in this repository as a submodule [found here](./kf-annotation-tools/workflows/kfdrc-somatic-snv-annot-workflow.cwl). To get the necessary files when downloading the repository, make sure use recursive flag for git clone: `git clone --recursive`. For details on the workflow and its components, please refer to the [workflow documentation](./kf-annotation-tools/docs/SOMATIC_SNV_ANNOT_README.md). The source repository for the annotation workflows can be [found here](https://github.com/kids-first/kf-annotation-tools).

  ### SV Callers

  - [Manta](https://github.com/Illumina/manta/tree/v1.4.0) `v1.4.0`, from Illumina, is used to call SVs. Output is also in VCF format, with calls filtered on `PASS`.

  #### SV Annotation

  - [AnnotSV](https://github.com/lgmgeo/AnnotSV/releases/tag/v3.1.1) `v3.1.1`, from Université de Strasbourg, is used to annotate the calls in the Manta PASS VCF.

  ### CNV Estimators

  - [Control-FREEC](https://github.com/BoevaLab/FREEC/tree/v11.6) `v11.6`, from the Bioinformatics Laboratory of Institut Curie and Boeva Lab, detects copy-number changes and allelic imbalances (including LOH)
     - The tool portion of the workflow is a port from the [Velsera](https://velsera.com/) team, with a slight tweak in image outputs
     - Also, the workflow wrapper limits what inputs and outputs are used based on our judgement of utility
  - [CNVkit](https://cnvkit.readthedocs.io/en/v0.9.3/) `v0.9.3`, from UCSF, infers and visualizes copy number
  - [THeTa2](https://github.com/kids-first/THetA/tree/v0.7.1) `v0.7.1`, from Brown University and CHOP D3b, is used to inform and adjust copy number calls from CNVkit with purity estimations
  - [GATK CNV](https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants) `v4.2.4.1`, from the Broad Institute, calls somatic CNVs using a Panel of Normals (PON) created using [this workflow](./workflows/kfdrc_gatk_create_cnv_pon_wf.cwl)
     - Read more about our the workflow and our PON recommendations in the [subworkflow doc](./docs/kfdrc_gatk_create_cnv_pon_wf.md)

  #### Copy Number Genotype Estimation Using B-Allele Frequency

  Control-FREEC and CNVkit can take advantage of B-Allele Frequency (BAF) integration for copy number genotype estimation and increased CNV accuracy. That is done using the following process:
  1. create a gVCF using either:
     - [Sentieon Alignment and gVCF Workflow](https://github.com/kids-first/kf-alignment-workflow/blob/master/workflows/kfdrc_sentieon_alignment_wf.cwl)
        - For details, please refer to the [workflow documentation](https://github.com/kids-first/kf-alignment-workflow/blob/master/docs/KFDRC_SENTIEON_ALIGNMENT_GVCF_WORKFLOW_README.md)
        - To run on CAVATICA, see our [public app](https://cavatica.sbgenomics.com/public/apps/cavatica/apps-publisher/kfdrc-sentieon-alignment-workflow)
     - [BWA Alignment and GATK HaplotypeCaller Workflow](https://github.com/kids-first/kf-alignment-workflow/blob/master/workflows/kfdrc_alignment_wf.cwl)
        - For details, please refer to the [workflow documentation](https://github.com/kids-first/kf-alignment-workflow/blob/master/readme.md)
        - To run on CAVATICA, see our [public app](https://cavatica.sbgenomics.com/public/apps/cavatica/apps-publisher/kfdrc-alignment-workflow)
  2. create a VCF by perform genotype refinement on the gVCF:
     - [Single Sample Genotyping Workflow](https://github.com/kids-first/kf-germline-workflow/blob/master/workflows/kfdrc-single-sample-genotyping-wf.cwl)
        - To run on CAVATICA, see our [public app](https://cavatica.sbgenomics.com/public/apps#cavatica/apps-publisher/kfdrc-single-sample-genotyping-wf/)
  3. Pass that VCF to the `b_allele` workflow input

  ### ecDNA Prediction
  - [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/documentation/GUIDE.md), from UCSD, is used to predict ecDNAs
     - Workflow documentation available [here](./docs/kfdrc_amplicon_architect.md)


  ## Running Whole Genome Sequencing (WGS) or Whole Exome Sequencing (WXS)

  The somatic variant workflow is designed to be able to process either WGS or
  WXS inputs. This functionality comes from usage of the `wgs_or_wxs` input enum.
  Depending on what is provided for this input, the tool will set the appropriate
  default values and check that the user has provided the correct inputs. For
  example, if the user sets the input to WGS the `lancet_padding` value will be
  set to 300; alternatively, if the user sets the input to WXS the
  `lancet_padding` value will be defaulted to 0. In either case, the user can
  override the defaults simply by providing their own value for `lancet_padding`
  in the inputs.

  The `wgs_or_wxs` flag also plays an integral role in the preparation of
  intervals throughout the pipeline.

  ## Calling Intervals, in short

  In order to perform variant calling, the user must provide some of these intervals:
  - `calling_regions`: genomic regions over which the SNV and SV callers will be run
  - `blacklist_regions`: genomic regions to remove from the SNV and SV calling regions
  - `coding_sequence_regions`: genomic regions covering the coding sequence regions of the genome to which the input BAM/CRAM/SAM was aligned
  - `cnv_blacklist_regions`: genomic regions to remove from the CNV calling regions

  The calling is the base calling regions for the sample. These regions are
  further modified by the two blacklist files. The blacklist regions first
  removed from the calling regions and the resulting regions are what are passed
  to callers. SNV and SV callers receive the regions that are `calling_regions`
  with `blacklist_regions` removed. CNV callers receive the regions that are
  `calling_regions` with `cnv_blacklist_regions` removed.

  In general, what is expected of the user to provide for these intervals is:
  - WGS
     - `calling_regions`: Whole genome calling regions (for GATK these are the regions of the FASTA on `chr1-22,X,Y` where there is no span of more than 200 Ns; for completeness we add `chrM`)
     - `blacklist_regions`: Regions that the user is confident are of no interest to SNV and SV calling, if any. Telomeres, centromeres, high signal regions, repetitive regions, etc.
     - `cnv_blacklist_regions`: Regions of the `calling_regions` the user is confident are of no interest to CNV calling, if any. These regions include chrM, telomeres, centromeres, high signal regions, repetitive regions, etc.
     - `coding_sequence_regions`: Regions of the genome that GENCODE defines as CDS; can be determined from [their GTF](https://www.gencodegenes.org/human/)
  - WXS
     - `calling_regions`: The experimental bait capture regions
     - `blacklist_regions`: Regions of the `calling_regions` that the user is confident are of no interest to SNV and SV calling, if any
     - `cnv_blacklist_regions` Regions of the capture region that the user is confident are of no interest to CNV calling, if any
     - `coding_sequence_regions`: This input is not used in WXS runs. Providing it does nothing

  The coding sequence regions are only relevant for calling WGS data using
  Lancet. Presenting Lancet with intervals spanning the whole genome will lead to
  extreme runtime and cost. By starting from an limited calling region and
  supplementing with calls from Mutect2 and/or Strelka2, Lancet is able to run on
  the scale of hours rather than days.

  There's a lot more to talk about with intervals, how they are prepared, and
  some finer details. If you are interested, see our [interval preparation walkthrough](./docs/interval_preparation_walkthrough.md).

  ## Partial Workflow Runs

  By default, the workflow will attempt to run every caller. However, the user can modify this behavior using various inputs.
  These inputs, and their behavior, are as follows:
  - `wgs_or_wxs`: WGS or WXS mode.
      - WXS mode disables AmpliconArchitect, even if `run_amplicon_architect` is set to true
      - Activates appropriate interval processing
  - `run_vardict`: Set to false to disable VarDict.
      - **Warning!** Enabling THeTa2 with VarDict disabled will throw an error!
  - `run_mutect2`: Set to false to disable Mutect2.
      - **Warning!** Enabling Lancet in WGS mode with Mutect2 and Strelka2 disabled will throw an error!
  - `run_strelka2`: Set to false to disable Strelka2.
      - **Warning!** Enabling Lancet in WGS mode with Strelka2 and Mutect2 disabled will throw an error!
  - `run_lancet`: Set to false to disable Lancet.
  - `run_controlfreec`: Set to false to disable Control-FREEC.
  - `run_cnvkit`: Set to false to disable CNVkit.
      - **Warning!** Enabling AmpliconArchitect with CNVkit disabled will throw an error!
      - **Warning!** Enabling THeTa2 with CNVkit disabled will throw an error!
  - `run_amplicon_architect`: Set to false to disable AmpliconArchitect.
  - `run_theta2`: Set to false to disable THeTa2.
  - `run_manta`: Set to false to disable Manta.
  - `run_gatk_cnv`: Set to false to disable GATK CNV.

  The first step of the workflow will check the user inputs and throw errors for impossible scenarios:
  - Enabling Lancet in WGS without enabling either Mutect2 or Strelka2 or user providing existing calls to the `lancet_input_vcf` input
  - Enabling AmpliconArchitect without enabling CNVkit or providing the Mosek config file
  - Enabling THeTa2 without enabling both VarDict and CNVkit
  - Enabling GATK CNV without providing a Panel of Normals

  # Tips to Run:

  1. For input cram files, be sure to have indexed them beforehand

  1. When in doubt, all of our reference files can be obtained from here: https://cavatica.sbgenomics.com/u/kfdrc-harmonization/kf-references/

  1. For Control-FREEC, it is highly recommended that you supply a VCF file with germline calls, GATK HaplotypeCaller recommended. Please also make sure the index for this file is available. Also, a range of input ploidy possibilities for the inputs are needed. You can simply use `2`, or put in a range, as an array, like 2, 3, 4. For mate orientation, you will need to specify, the drop down and tool doc explains your options.

  1. As a CAVATICA app, default references for hg38 are already pre-populated, as well as some default settings - i.e., number of threads, coefficient of variation input for Control-FREEC, and `PASS` filter tool mode.

  1. If running Lancet for WGS, and you have existing variant calls that you'd like to use from a previous Strelka2 and/or Mutect2 run for example, provide them to `lancet_input_vcf` input. Otherwise you _must_ also enable `run_strelka2` and/or `run_mutect2`

  1. `select_vars_mode`: On occasion, using GATK's `SelectVariants` tool will fail, so a simple `grep` mode on `PASS` can be used instead.
  Related, `bcftools_filter_vcf` is built in as a convenience in case your b allele frequency file has not been filtered on `PASS`.
  You can use the `include_expression` `Filter="PASS"` to achieve this.

  1. The `SM:` sample IDs in the input alignment files **must** match `input_tumor_name` and `input_normal_name`. If you need the output VCF sample names to differ from the input, enter what is currently in `SM:` into the `old_tumor_name` and `old_normal_name` fields, and enter the desired sample names in `input_tumor_name` and `input_normal_name`. The workflow will rename the samples in the VCF for you.

  1. Our reference inputs can be obtained with an account [from CAVATICA](https://cavatica.sbgenomics.com/u/kfdrc-harmonization/kf-references/). Below are some directions on how to make them:
     - `reference_fasta`: [Homo_sapiens_assembly38.fasta](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK
     - `reference_dict`: [Homo_sapiens_assembly38.dict](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK
     - `calling_regions`: [wgs_calling_regions.hg38.interval_list](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK. **To create our canonical calling intervals, edit this file by leaving only entries related to chr1-22,X,Y,M. M may need to be added.**
     - `cnv_blacklist_regions`: `somatic-hg38_CNV_and_centromere_blacklist.hg38liftover.list` Blacklist regions that include centromeres to exclude from CNV calling
     - `coding_sequence_regions`: `GRCh38.gencode.v31.CDS.merged.bed` For Lancet WGS, it's highly recommended to use CDS bed as the starting point and supplement with the regions of calls from Strelka2 & Mutect2. Our CDS regions were obtained from GENCODE, [release 31](https://www.gencodegenes.org/human/release_31.html) using this GTF file [gencode.v31.primary_assembly.annotation.gtf.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.primary_assembly.annotation.gtf.gz) and parsing features for `UTR`, `start codon`, `stop codon`, and `exon`, then using bedtools sort and merge after converting coordinates into bed format.
     - `cnvkit_annotation_file`: [refFlat_HG38.txt](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz) gunzip this file from UCSC
     - `af_only_gnomad_vcf`: [af-only-gnomad.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38) - need a valid google account, this is a link to the best practices google bucket from Broad GATK.
     - `exac_common_vcf`: [small_exac_common_3.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38) - need a valid google account, this is a link to the best practices google bucket from Broad GATK.
     - `vep_cache`: `homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz` from [VEP release 105](https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/) - variant effect predictor cache. Current production workflow uses this version.
     - `genomic_hotspots`: `tert.bed`. Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots
     - `protein_snv_hotspots`: [kfdrc_protein_snv_cancer_hotspots_20240718.txt](https://www.cancerhotspots.org/files/hotspots_v2.xls). Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots. Recommend pulling the two relevant columns for SNVs only, and convert to TSV
     - `protein_indel_hotspots`: [protein_indel_cancer_hotspots_v2.ENS105_liftover.tsv](https://www.cancerhotspots.org/files/hotspots_v2.xls). Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots.
     - `echtvar_anno_zips`: `gnomad_3_1_1.vwb_subset.echtvar_0_1_9.zip`.Echtvar annotation ZIP file. For more info on how this file was created, see our [Ectvar documentation](./docs/echtvar_gnomad_annotation.md).
     - `custom_enst`: `kf_isoform_override.tsv`. As of VEP 104, several genes have had their canonical transcripts redefined. While the VCF will have all possible isoforms, this affects MAF file output and may results in representative protein changes that defy historical expectations
     - `funcotator_data_sources_tgz`: [funcotator_dataSources.v1.6.20190124s.tar.gz](https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator) - need a valid google account, this is a link to pre-packaged datasources from the Broad Institute. Any of the tar.gz files will do.
     - `annotsv_annotations_dir_tgz`: [annotsv_311_annotations_dir.tgz] - These annotations are simply those from the install-human-annotation installation process run during AnnotSV installation. Specifically these are the annotations installed with v3.1.1 of the software. Newer or older annotations can be slotted in here as needed.
     - `aa_data_repo`: [GRCh38_indexed.tar.gz](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/) - AmpliconArchitect reference TAR ball. Only needed/used if input is WGS
     - `mosek_license_file`: [mosek.lic](https://www.mosek.com/license/request/) - required if AmpliconArchitect is run. License is good for one year and is renewable

  1. There are some flags that we set as defaults. Feel free to alter these if your experiment calls for something else:
     - `cfree_coeff_var`: 0.05
     - `cfree_contamination_adjustment`: FALSE
     - `bcftools_public_filter`: `'FILTER="PASS"|INFO/HotSpotAllele=1'`. This phrase will allow `PASS` only **or** `HotSpotAllele` variants into the public version of variant call output.
     - `disable_hotspot_annotation`: false
     - `maf_center`: `"."`. Sequencing center of variant called
     - `aa_data_ref_version`: `"GRCh38"`. Genome reference version used

  1. There are some flags that need to be set by the user at runtime. These are our recommendations:
     - `gatk_filter_name`: `["NORM_DP_LOW", "GNOMAD_AF_HIGH"]`. These correspond to the recommended filter expression.
     - `gatk_filter_expression`: `["vc.getGenotype('<input_normal_name> ').getDP() <= 7"), "gnomad_3_1_1_AF != '.' && gnomad_3_1_1_AF > 0.001  && gnomad_3_1_1_FILTER=='PASS'"]`. Array of filter expressions to establish criteria to tag variants with. See [annotation subworkflow docs](./docs/kfdrc_annotation_subworkflow.md) for a more detailed explanation. For more information on filter expressions, see the [GATK JEXL docs](https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration).
     - `extra_arg`: This can be used to add special params to Strelka2. It is currently more of an "unsticking param". For edge cases where strelka2 seems to hang, setting this to `--max-input-depth 10000` can balance performance and consistency in results

  1. Output files (Note, all VCF files that don't have an explicit index output have index files output as as secondary file.  In other words, they will be captured at the end of the workflow):
      - Simple variant callers
          - Strelka2
              - `strelka2_prepass_vcf`: Combined SNV + INDEL file with renamed Sample IDs. Has all soft `FILTER` values generated by variant caller
              - `strelka2_protected_outputs`: Array of files containing MAF format of PASS hits, `PASS` VCF with annotation pipeline soft `FILTER`-added values, and VCF index
              - `strelka2_public_outputs`: Same as above, except MAF and VCF have had entries with soft `FILTER` values removed
          - Mutect2
              - `mutect2_filtered_vcf`: VCF with SNV, MNV, and INDEL variant calls. Contains all soft `FILTER` values generated by variant caller. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter
              - `mutect2_protected_outputs`: Array of files containing MAF format of PASS hits, `PASS` VCF with annotation pipeline soft `FILTER`-added values, and VCF index
              - `mutect2_public_outputs`: Same as above, except MAF and VCF have had entries with soft `FILTER` values removed
          - VarDict Java
              - `vardict_prepass_vcf`: VCF with SNV, MNV, and INDEL variant calls. Contains all soft `FILTER` values generated by variant caller. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter
              - `vardict_protected_outputs`: Array of files containing MAF format of PASS hits, `PASS` VCF with annotation pipeline soft `FILTER`-added values, and VCF index
              - `vardict_public_outputs`: Same as above, except MAF and VCF have had entries with soft `FILTER` values removed
          - Lancet
              - `lancet_prepass_vcf`: VCF with SNV, MNV, and INDEL variant calls. Contains all soft `FILTER` values generated by variant caller, Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter
              - `lancet_protected_outputs`: Array of files containing MAF format of PASS hits, `PASS` VCF with annotation pipeline soft `FILTER`-added values, and VCF index
              - `lancet_public_outputs`: Same as above, except MAF and VCF have had entries with soft `FILTER` values removed
      - Structural variant callers
          - Manta
              - `manta_pass_vcf`: SV call filtered on `PASS`, from manta
              - `manta_pass_tbi`: Index file of above bgzipped VCF
              - `manta_prepass_vcf`: SV results with all `FILTER` categories for manta. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter
              - `annotsv_annotated_calls`: This file contains all records from the `manta_pass_vcf` that AnnotSV could annotate.
              - `annotsv_unannotated_calls`: This file contains all records from the `manta_pass_vcf` that AnnotSV could not annotate.
      - Copy number variation callers
          - ControlFREEC
              - `ctrlfreec_pval`: CNV calls with copy number and p value confidence, a qualitative "gain, loss, neutral" assignment, and genotype with uncertainty assigned from Control-FREEC.  See author manual for more details.
              - `ctrlfreec_config`: Config file used to run Control-FREEC.  Has some useful information on what parameters were used to run the tool.
              - `ctrlfreec_pngs`: Plots of BAF log2 ratio and ratio of tumor/normal copy number coverage.  Pink line in the middle of ratio plots is the median ratio.
              - `ctrlfreec_bam_ratio`: BAM ratio text file.  Contain ratio, median ratio (used to inform `ctrlfreec_pval`), CNV estimates, BAF estimate, and genotype estimate.
              - `ctrlfreec_bam_seg`: In-house generated SEG file based on ratio file.  Provided as a convenience for compatibility with tools that require inputs in legacy microarray format.
              - `ctrlfreec_baf`: BAF estimations.
              - `ctrlfreec_info`: Contains useful run time information, like ploidy used for analysis, and window size
          - CNVkit
              - `cnvkit_cnr`: Copy number ratio
              - `cnvkit_cnn_output`: Normal/control sample copy number
              - `cnvkit_calls`: Tumor/sample copy number
              - `cnvkit_metrics`: Basic SEG count and call stats
              - `cnvkit_gainloss`: Per-gene log2 ratio
              - `cnvkit_seg`: Classic microarray-style SEG file
              - `theta2_calls`: Purity-adjusted CNVkit copy number calls based on THeTa2 results
              - `theta2_seg`: Purity-adjusted CNVkit SEG file based on THeTa2 results
              - `theta2_subclonal_results`: THeTa2 subclone purity results
              - `theta2_subclonal_cns`: THeTa2 sublone CNS
              - `theta2_subclone_seg`: THeTa2 subclone SEG file
          - GATK CNV
              - `gatk_copy_ratio_segments_tumor`: Called copy-ratio-segments file for the tumor sample. This is a tab-separated values (TSV) file with SAM-style header containing a read group sample name, a sequence dictionary, a row specifying the column headers output by `CallCopyRatioSegments` and the corresponding entry rows.
              - `gatk_copy_ratio_segments_normal`: Called copy-ratio-segments file for the normal sample. This is a tab-separated values (TSV) file with SAM-style header containing a read group sample name, a sequence dictionary, a row specifying the column headers output by `CallCopyRatioSegments` and the corresponding entry rows.
              - `gatk_cnv_denoised_tumor_plot`: Denoised-plot file that covers the entire range of the copy ratios
              - `gatk_cnv_denoised_normal_plot`: Denoised-plot file that covers the entire range of the copy ratios
              - `gatk_cnv_funcotated_called_file_tumor`: TSV where each row is a segment and the annotations are the covered genes and which genes+exon is overlapped by the segment breakpoints.
              - `gatk_cnv_funcotated_called_gene_list_file_tumor`: TSV where each row is a gene and the annotations are the covered genes and which genes+exon is overlapped by the segment breakpoints.
      - ecDNA Prediction
          - `aa_summary`: Summary for all amplicons detected by AmpliconArchitect
          - `aa_cycles`: Text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant source) and their copy counts
          - `aa_graph`: A text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
          - `aa_sv_png`: PNG image file displaying the SV view of AmpliconArchitect
          - `aa_classification_profiles`: abstract classification of the amplicon

  1. Docker images - the workflow tools will automatically pull them, but as a convenience are listed below:
      - `Strelka2`: pgc-images.sbgenomics.com/d3b-bixu/strelka
      - `add common fields to Strelka2`: pgc-images.sbgenomics.com/d3b-bixu/add-strelka2-fields:1.0.0
      - `Mutect2` and `GATK SNV` tools: pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
      - `GATK CNV`: broadinstitute/gatk:4.2.4.1
      - `Lancet`: pgc-images.sbgenomics.com/d3b-bixu/lancet:1.0.7
      - `VarDict Java`: pgc-images.sbgenomics.com/d3b-bixu/vardict:1.7.0
      - `Control-FREEC`: images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1
      - `CNVkit`: images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3
      - `THetA2`: pgc-images.sbgenomics.com/d3b-bixu/theta2:0.7.1
      - `samtools`: pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9
      - `Variant Effect Predictor`: ensemblorg/ensembl-vep:release_105.0
      - `Kids Fist VCF2MAF`: pgc-images.sbgenomics.com/d3b-bixu/kf_vcf2maf:v1.0.3
      - `Manta`: pgc-images.sbgenomics.com/d3b-bixu/manta:1.4.0
      - `bcftools` and `vcftools`: pgc-images.sbgenomics.com/d3b-bixu/bvcftools:latest
      - `annotsv`: pgc-images.sbgenomics.com/d3b-bixu/annotsv:3.1.1
      - `AmpliconArchitect`: jluebeck/prepareaa:v0.1203.10
      - `Echtvar`: pgc-images.sbgenomics.com/d3b-bixu/echtvar:0.1.9

  1. For highly complex samples, some tools have shown themselves to require memory allocation adjustments:
     Manta, GATK LearnReadOrientationModel, GATK GetPileupSummaries, GATK FilterMutectCalls and VarDict.
     Optional inputs exist to expand the memory allocation for these jobs: manta_memory, learnorientation_memory,
     getpileup_memory, filtermutectcalls_memory, and vardict_ram, respectively.
     For the java tools (VarDict and GATK), these values represent limits on the memory that can
     be used for the respective jobs. Tasks will go to these values and not exceed it. They may succeed or fail,
     but they will not exceed the limit established. The memory allocations for these is hardcapped. The memory
     allocation option for Manta, conversely, is a soft cap. The memory requested will be allocated for the job
     on a particular machine but once the task is running on the machine it may exceed that requested value. For example,
     if Manta's memory allocation is set to 10 GB it will have 10 GB allocated to it at task creation, but, if the
     task ends up running on a machine with more memory available, the task may use it. Setting a value here for Manta
     will not prevent Manta from taking more than that value. The memory usage in Manta is limited by the machine hardware.
     As such the option for Manta memory allocation is described as soft cap. For more information on Manta resource
     usage see their [documentation](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#runtime-hardware-requirements).

  1. The optional `b_allele` file can be generated using our [Single Sample Genotyping Workflow](#copy-number-genotype-estimation-using-b-allele-frequency).

  ## Other Resources
  - dockerfiles: https://github.com/d3b-center/bixtools
requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
inputs:
  # Required
  indexed_reference_fasta: {type: 'File', secondaryFiles: [{pattern: ".fai", required: true}, {pattern: "^.dict", required: true}],
    "sbg:suggestedValue": {class: File, path: 60639014357c3a53540ca7a3, name: Homo_sapiens_assembly38.fasta, secondaryFiles: [{class: File,
          path: 60639016357c3a53540ca7af, name: Homo_sapiens_assembly38.fasta.fai}, {class: File, path: 60639019357c3a53540ca7e7,
          name: Homo_sapiens_assembly38.dict}]}}
  input_tumor_aligned:
    type: File
    secondaryFiles: [{pattern: ".bai", required: false}, {pattern: "^.bai", required: false}, {pattern: ".crai", required: false},
      {pattern: "^.crai", required: false}]
    doc: "tumor BAM or CRAM"
  input_tumor_name: {type: string, doc: "Desired sample name for tumor in output VCFs"}
  old_tumor_name: {type: 'string?', doc: "If `SM:` sample name in te align file is different than `input_tumor_name`, you **must**
      provide it here"}
  input_normal_aligned:
    type: File
    secondaryFiles: [{pattern: ".bai", required: false}, {pattern: "^.bai", required: false}, {pattern: ".crai", required: false},
      {pattern: "^.crai", required: false}]
    doc: "normal BAM or CRAM"
  input_normal_name: {type: string, doc: "Desired sample name for normal in output VCFs"}
  old_normal_name: {type: 'string?', doc: "If `SM:` sample name in te align file is different than `input_normal_name`, you **must**
      provide it here"}
  calling_regions: {type: 'File', doc: "BED or INTERVALLIST file containing a set of genomic regions over which the callers will be
      run. For WGS, this should be the wgs_calling_regions.interval_list. For WXS, the user must provide the appropriate regions for
      their analysis."}
  blacklist_regions: {type: 'File?', doc: "BED or INTERVALLIST file containing a set of genomic regions to remove from the calling
      regions for SNV and SV calling."}
  cnv_blacklist_regions: {type: 'File?', doc: "BED or INTERVALLIST file containing a set of genomic regions to remove from the calling
      regions for CNV calling only!", "sbg:suggestedValue": {class: File, path: 663d2bcc27374715fccd8c6d, name: somatic-hg38_CNV_and_centromere_blacklist.hg38liftover.list}}
  coding_sequence_regions: {type: 'File?', doc: "BED or INTERVALLIST file containing the coding sequence regions for the provided
      reference. This input is used to create custom intervals for WGS Lancet Calling.", "sbg:suggestedValue": {class: File, path: 5f500135e4b0370371c051c0,
      name: GRCh38.gencode.v31.CDS.merged.bed}}
  cfree_ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try"}
  cnvkit_annotation_file: {type: 'File', doc: "refFlat.txt file", "sbg:suggestedValue": {class: File, path: 5f500135e4b0370371c051c1,
      name: refFlat_HG38.txt}}
  extra_arg: {type: 'string?', doc: "Add special options to config file, i.e. --max-input-depth 10000"}
  strelka2_cores: {type: 'int?', doc: "Adjust number of cores used to run strelka2", default: 18}
  mutect2_af_only_gnomad_vcf: {type: 'File', secondaryFiles: [{pattern: ".tbi", required: true}], "sbg:suggestedValue": {class: File,
      path: 5f50018fe4b054958bc8d2e3, name: af-only-gnomad.hg38.vcf.gz, secondaryFiles: [{class: File, path: 5f50018fe4b054958bc8d2e5,
          name: af-only-gnomad.hg38.vcf.gz.tbi}]}}
  mutect2_exac_common_vcf: {type: 'File', secondaryFiles: [{pattern: ".tbi", required: true}], "sbg:suggestedValue": {class: File,
      path: 5f500135e4b0370371c051ad, name: small_exac_common_3.hg38.vcf.gz, secondaryFiles: [{class: File, path: 5f500135e4b0370371c051af,
          name: small_exac_common_3.hg38.vcf.gz.tbi}]}}
  output_basename: {type: 'string', doc: "String value to use as basename for outputs"}
  wgs_or_wxs: {type: {type: enum, name: wgs_or_wxs, symbols: ["WGS", "WXS"]}, doc: "Select if this run is WGS or WXS"}
  count_panel_of_normals: {type: 'File?', doc: "Path to read-count PoN created by the panel workflow. Significantly reduces quality
      of calling if not provided!", "sbg:fileTypes": "HDF5"}
  run_funcotatesegments: {type: 'boolean?', default: true, doc: "If true, run Funcotator on the called copy-ratio segments. This will
      generate both a simple TSV and a gene list."}
  funcotator_data_sources_tgz: {type: 'File?', doc: "Path to tar.gz containing the data sources for Funcotator to create annotations.",
    "sbg:fileTypes": "TAR, TAR.GZ, TGZ", "sbg:suggestedValue": {class: File, path: 60e5f8636a504e4e0c6408d8, name: funcotator_dataSources.v1.6.20190124s.tar.gz}}
  funcotator_minimum_segment_size: {type: 'int?', doc: "The minimum number of bases for a variant to be annotated as a segment. Recommended
      to be changed only for use with FuncotateSegments. If you encounter 'Variant context does not represent a copy number segment'
      error, set this value lower than the length of the failed segment."}
  aa_data_repo: {type: 'File?', doc: "Reference tar ball obtained from https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/",
    "sbg:suggestedValue": {class: File, path: 62fcf4d40d34597148589e14, name: GRCh38_indexed.tar.gz}}
  aa_data_ref_version: {type: ['null', {type: enum, name: aa_data_ref_version, symbols: ["GRCh38", "hg19", "GRCh37", "mm10", "GRCm38"]}],
    doc: "Genome version in data repo to use", default: "GRCh38"}
  mosek_license_file: {type: 'File?', doc: "This tool uses some software that requires a license file. Only provide if input is WGS.
      You can get a personal or institutional one from https://www.mosek.com/license/request/.", "sbg:suggestedValue": {class: File,
      path: 62fcf4d40d34597148589e15, name: mosek.lic}}
  run_vardict: {type: 'boolean?', default: true, doc: "Set to false to disable Vardict. Warning: Vardict is required to run Theta2!"}
  run_mutect2: {type: 'boolean?', default: true, doc: "Set to false to disable Mutect2. Warning: Mutect2 is required to run Lancet
      in WGS mode!"}
  run_strelka2: {type: 'boolean?', default: true, doc: "Set to false to disable Strelka2. Warning: Strelka2 is required to run Lancet
      in WGS mode!"}
  run_lancet: {type: 'boolean?', default: true, doc: "Set to false to disable Lancet."}
  lancet_input_vcf: {type: 'File[]?', doc: "Use only if you need to re-run lancet and have pre-existing strelka2 and/or mutect2 results"}
  run_controlfreec: {type: 'boolean?', default: true, doc: "Set to false to disable ControlFreeC."}
  run_cnvkit: {type: 'boolean?', default: true, doc: "Set to false to disable CNVkit. Warning: CNVkit is required to run both Amplicon
      Architect and Theta2!"}
  run_amplicon_architect: {type: 'boolean?', default: true, doc: "Set to false to disable Amplicon Architect."}
  run_theta2: {type: 'boolean?', default: true, doc: "Set to false to disable Theta2."}
  run_manta: {type: 'boolean?', default: true, doc: "Set to false to disable Manta."}
  run_gatk_cnv: {type: 'boolean?', default: true, doc: "Set to false to disable GATK CNV."}
  run_calmd_bam: { type: 'boolean?', default: false, doc: "Override logic to skip calmd when input is not cram" }
  annotsv_annotations_dir_tgz: {type: 'File?', doc: "TAR.GZ'd Directory containing annotations for AnnotSV", "sbg:fileTypes": "TAR,
      TAR.GZ, TGZ", "sbg:suggestedValue": {class: File, path: 6328ab26d01163633dabcc2e, name: annotsv_311_plus_ens105_annotations_dir.tgz}}
  cfree_threads: {type: 'int?', default: 16, doc: "For ControlFreeC. Recommend 16 max, as I/O gets saturated after that losing any
      advantage"}
  cfree_mate_orientation_control: {type: ['null', {type: enum, name: mate_orientation_control, symbols: ["0", "FR", "RF", "FF"]}],
    default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  cfree_mate_orientation_sample: {type: ['null', {type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"]}], default: "FR",
    doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  lancet_ram: {type: 'int?', default: 12, doc: "Adjust in rare circumstances in which 12 GB is not enough"}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], default: "gatk", doc: "Choose
      'gatk' for SelectVariants tool, or 'grep' for grep expression"}
  min_theta2_frac: {type: 'float?', default: 0.01, doc: "Minimum fraction of genome with copy umber alterations.  Default is 0.05,
      recommend 0.01"}
  vardict_cpus: {type: 'int?', default: 8, doc: "Number of CPUs for Vardict to use"}
  vardict_min_vaf: {type: 'float?', default: 0.05, doc: "Min variant allele frequency for vardict to consider. Recommend 0.05"}
  vardict_ram: {type: 'int?', default: 16, doc: "GB of RAM to allocate to Vardict (hard-capped)"}
  exome_flag: {type: 'string?', doc: "Whether to run in exome mode for callers. Y for WXS, N for WGS"}
  lancet_window: {type: 'int?', doc: "Window size for lancet.  Recommend 500 for WGS; 600 for exome+"}
  lancet_padding: {type: 'int?', doc: "Recommend 0 if interval file padded already, half window size if not. Recommended: 0 for WXS;
      300 for WGS"}
  vardict_padding: {type: 'int?', doc: "Padding to add to input intervals, recommend 0 if intervals already padded such as in WXS,
      150 if not such as in WGS"}
  cnvkit_wgs_mode: {type: 'string?', doc: "for WGS mode, input Y. leave blank for WXS/hybrid mode"}
  i_flag: {type: 'string?', doc: "Flag to intersect germline calls on padded regions. Use N if you want to skip this or have a WGS
      run"}
  b_allele: {type: 'File?', secondaryFiles: [{pattern: ".tbi", required: true}], doc: "germline calls, needed for BAF.  GATK HC VQSR
      input recommended.  Tool will prefilter for germline and pass if expression given"}
  cfree_coeff_var: {type: 'float?', default: 0.05, doc: "Coefficient of variation to set window size.  Default 0.05 recommended"}
  cfree_contamination_adjustment: {type: 'boolean?', doc: "TRUE or FALSE to have ControlFreec estimate normal contam"}
  cfree_sex: {type: ['null', {type: enum, name: cfree_sex, symbols: ["XX", "XY"]}], doc: "If known, XX for female, XY for male", default: "XX"}
  cnvkit_sex: {type: ['null', {type: enum, name: cnvkit_sex, symbols: ["x", "y"]}], doc: "Sex, for simplicity x for female y for male",
    default: "x"}
  combined_include_expression: {type: 'string?', doc: "Theta2 Purity value: Filter expression if vcf has non-PASS combined calls,
      use as-needed, default set for VarDict Java for VarDict", default: FILTER="PASS" && (INFO/STATUS="Germline" | INFO/STATUS="StrongSomatic")}
  combined_exclude_expression: {type: 'string?', doc: "Theta2 Purity value: Filter expression if vcf has non-PASS combined calls,
      use as-needed"}
  use_manta_small_indels: {type: 'boolean?', default: false, doc: "Should the program use the small indels output from Manta in Strelka2
      calling?"}
  learnorientation_memory: {type: 'int?', doc: "GB of memory to allocate to GATK LearnReadOrientationModel; defaults to 4 (hard-capped)"}
  getpileup_memory: {type: 'int?', doc: "GB of memory to allocate to GATK GetPileupSummaries; defaults to 2 (hard-capped)"}
  filtermutectcalls_memory: {type: 'int?', doc: "GB of memory to allocate to GATK FilterMutectCalls; defaults to 4 (hard-capped)"}
  manta_memory: {type: 'int?', doc: "GB of memory to allocate to Manta; defaults to 10 (soft-capped)"}
  manta_cores: {type: 'int?', doc: "Number of cores to allocate to Manta; defaults to 18"}
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache", "sbg:suggestedValue": {class: File, path: 6332f8e47535110eb79c794f,
      name: homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz}}
  vep_ram: {type: 'int?', default: 32, doc: "In GB, may need to increase this value depending on the size/complexity of input"}
  vep_cores: {type: 'int?', default: 16, doc: "Number of cores to use. May need to increase for really large inputs"}
  vep_buffer_size: {type: 'int?', default: 5000, doc: "Increase or decrease to balance speed and memory usage"}
  dbnsfp: {type: 'File?', secondaryFiles: [.tbi, ^.readme.txt], doc: "VEP-formatted plugin file, index, and readme file containing
      dbNSFP annotations"}
  dbnsfp_fields: {type: 'string?', doc: "csv string with desired fields to annotate. Use ALL to grab all", default: 'clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,Interpro_domain,GTEx_V8_gene,GTEx_V8_tissue'}
  merged: {type: 'boolean?', doc: "Set to true if merged cache used", default: true}
  cadd_indels: {type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD indel annotations"}
  cadd_snvs: {type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD SNV annotations"}
  run_cache_existing: {type: 'boolean?', doc: "Run the check_existing flag for cache"}
  run_cache_af: {type: 'boolean?', doc: "Run the allele frequency flags for cache"}
  strelka2_retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep, i.e. for strelka2 `MQ,MQ0,QSI,HotSpotAllele`",
    default: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,gnomad_3_1_1_AF_non_cancer_afr,gnomad_3_1_1_AF_non_cancer_ami,gnomad_3_1_1_AF_non_cancer_asj,gnomad_3_1_1_AF_non_cancer_eas,gnomad_3_1_1_AF_non_cancer_fin,gnomad_3_1_1_AF_non_cancer_mid,gnomad_3_1_1_AF_non_cancer_nfe,gnomad_3_1_1_AF_non_cancer_oth,gnomad_3_1_1_AF_non_cancer_raw,gnomad_3_1_1_AF_non_cancer_sas,gnomad_3_1_1_AF_non_cancer_amr,gnomad_3_1_1_AF_non_cancer_popmax,gnomad_3_1_1_AF_non_cancer_all_popmax,gnomad_3_1_1_FILTER,MQ,MQ0,QSI,HotSpotAllele"}
  strelka2_retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  strelka2_retain_ann: {type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF",
    default: "HGVSg"}
  mutect2_retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep, i.e. for mutect2 `MBQ,TLOD,HotSpotAllele`",
    default: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,gnomad_3_1_1_AF_non_cancer_afr,gnomad_3_1_1_AF_non_cancer_ami,gnomad_3_1_1_AF_non_cancer_asj,gnomad_3_1_1_AF_non_cancer_eas,gnomad_3_1_1_AF_non_cancer_fin,gnomad_3_1_1_AF_non_cancer_mid,gnomad_3_1_1_AF_non_cancer_nfe,gnomad_3_1_1_AF_non_cancer_oth,gnomad_3_1_1_AF_non_cancer_raw,gnomad_3_1_1_AF_non_cancer_sas,gnomad_3_1_1_AF_non_cancer_amr,gnomad_3_1_1_AF_non_cancer_popmax,gnomad_3_1_1_AF_non_cancer_all_popmax,gnomad_3_1_1_FILTER,MBQ,TLOD,HotSpotAllele"}
  mutect2_retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  mutect2_retain_ann: {type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF",
    default: "HGVSg"}
  lancet_retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep, i.e. for lancet `MS,FETS,HotSpotAllele`",
    default: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,gnomad_3_1_1_AF_non_cancer_afr,gnomad_3_1_1_AF_non_cancer_ami,gnomad_3_1_1_AF_non_cancer_asj,gnomad_3_1_1_AF_non_cancer_eas,gnomad_3_1_1_AF_non_cancer_fin,gnomad_3_1_1_AF_non_cancer_mid,gnomad_3_1_1_AF_non_cancer_nfe,gnomad_3_1_1_AF_non_cancer_oth,gnomad_3_1_1_AF_non_cancer_raw,gnomad_3_1_1_AF_non_cancer_sas,gnomad_3_1_1_AF_non_cancer_amr,gnomad_3_1_1_AF_non_cancer_popmax,gnomad_3_1_1_AF_non_cancer_all_popmax,gnomad_3_1_1_FILTER,MS,FETS,HotSpotAllele"}
  lancet_retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  lancet_retain_ann: {type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF",
    default: "HGVSg"}
  vardict_retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep, i.e. for consensus `MSI,MSILEN,SOR,SSF,HotSpotAllele`",
    default: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,gnomad_3_1_1_AF_non_cancer_afr,gnomad_3_1_1_AF_non_cancer_ami,gnomad_3_1_1_AF_non_cancer_asj,gnomad_3_1_1_AF_non_cancer_eas,gnomad_3_1_1_AF_non_cancer_fin,gnomad_3_1_1_AF_non_cancer_mid,gnomad_3_1_1_AF_non_cancer_nfe,gnomad_3_1_1_AF_non_cancer_oth,gnomad_3_1_1_AF_non_cancer_raw,gnomad_3_1_1_AF_non_cancer_sas,gnomad_3_1_1_AF_non_cancer_amr,gnomad_3_1_1_AF_non_cancer_popmax,gnomad_3_1_1_AF_non_cancer_all_popmax,gnomad_3_1_1_FILTER,MSI,MSILEN,SOR,SSF,HotSpotAllele"}
  vardict_retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  vardict_retain_ann: {type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF",
    default: "HGVSg"}
  genomic_hotspots: {type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to
      hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a423, name: tert.bed}]}
  protein_snv_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid
      positions corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 66980e845a58091951d53984, name: kfdrc_protein_snv_cancer_hotspots_20240718.txt}]}
  protein_indel_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino
      acid position ranges corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 663d2bcc27374715fccd8c6f, name: protein_indel_cancer_hotspots_v2.ENS105_liftover.tsv}]}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to
      create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  echtvar_anno_zips: {type: 'File[]?', doc: "Annotation ZIP files for echtvar anno", "sbg:suggestedValue": [{class: File, path: 65c64d847dab7758206248c6,
        name: gnomad.v3.1.1.custom.echtvar.zip}]}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"\
      ]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration,
      recommend: [`vc.getGenotype('inputs.input_normal_name').getDP() <= 7)`, `gnomad_3_1_1_AF > 0.001`]"}
  disable_hotspot_annotation: {type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task.", default: false}
  disable_vep_annotation: {type: 'boolean?', doc: "Disable VEP Annotation and skip this task.", default: false}
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}
  custom_enst: {type: 'File?', doc: "Use a file with ens tx IDs for each gene to override VEP PICK", "sbg:suggestedValue": {class: File,
      path: 663d2bcc27374715fccd8c65, name: kf_isoform_override.tsv}}
outputs:
  aa_summary: {type: 'File?', doc: "summary for all amplicons detected by AA", outputSource: amplicon_architect/aa_summary}
  aa_cycles: {type: 'File[]?', doc: "text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence,
      discordant, concordant, source) and their copy counts", outputSource: amplicon_architect/aa_cycles}
  aa_graph: {type: 'File[]?', doc: 'A text file for each amplicon listing the edges in the breakpoint graph, their categorization
      (sequence, discordant, concordant, source) and their copy counts', outputSource: amplicon_architect/aa_graph}
  aa_sv_png: {type: 'File[]?', doc: "PNG image file displaying the SV view of AA", outputSource: amplicon_architect/aa_sv_png}
  aa_classification_profiles: {type: 'File[]?', doc: "abstract classification of the amplicon", outputSource: amplicon_architect/aa_classification_profiles}
  aa_gene_list: {type: 'File[]?', doc: "genes present on amplicons with each classification", outputSource: amplicon_architect/aa_gene_list}
  ctrlfreec_pval: {type: 'File?', outputSource: controlfreec/ctrlfreec_pval}
  ctrlfreec_config: {type: 'File?', outputSource: controlfreec/ctrlfreec_config}
  ctrlfreec_pngs: {type: 'File[]?', outputSource: controlfreec/ctrlfreec_pngs}
  ctrlfreec_bam_ratio: {type: 'File?', outputSource: controlfreec/ctrlfreec_bam_ratio}
  ctrlfreec_bam_seg: {type: 'File?', outputSource: controlfreec/ctrlfreec_bam_seg}
  ctrlfreec_baf: {type: 'File?', outputSource: controlfreec/ctrlfreec_baf}
  ctrlfreec_info: {type: 'File?', outputSource: controlfreec/ctrlfreec_info}
  cnvkit_cnr: {type: 'File?', outputSource: cnvkit/cnvkit_cnr}
  cnvkit_cnn_output: {type: 'File?', outputSource: cnvkit/cnvkit_cnn_output}
  cnvkit_calls: {type: 'File?', outputSource: cnvkit/cnvkit_calls}
  cnvkit_metrics: {type: 'File?', outputSource: cnvkit/cnvkit_metrics}
  cnvkit_gainloss: {type: 'File?', outputSource: cnvkit/cnvkit_gainloss}
  cnvkit_seg: {type: 'File?', outputSource: cnvkit/cnvkit_seg}
  cnvkit_scatter_plot: {type: 'File?', outputSource: cnvkit/cnvkit_scatter_plot}
  cnvkit_diagram: {type: 'File?', outputSource: cnvkit/cnvkit_diagram}
  theta2_calls: {type: 'File?', outputSource: theta2_purity/theta2_adjusted_cns}
  theta2_seg: {type: 'File?', outputSource: theta2_purity/theta2_adjusted_seg}
  theta2_subclonal_results: {type: 'File[]?', outputSource: expression_flatten_subclonal_results/output}
  theta2_subclonal_cns: {type: 'File[]?', outputSource: theta2_purity/theta2_subclonal_cns}
  theta2_subclone_seg: {type: 'File[]?', outputSource: theta2_purity/theta2_subclone_seg}
  strelka2_public_outputs: {type: 'File[]', outputSource: strelka2/strelka2_public_outputs}
  strelka2_protected_outputs: {type: 'File[]', outputSource: strelka2/strelka2_protected_outputs}
  strelka2_prepass_vcf: {type: 'File', outputSource: strelka2/strelka2_prepass_vcf}
  manta_pass_vcf: {type: 'File?', outputSource: manta/manta_pass_vcf}
  manta_prepass_vcf: {type: 'File?', outputSource: manta/manta_prepass_vcf}
  annotsv_annotated_calls: {type: 'File?', outputSource: manta/annotsv_annotated_calls}
  annotsv_unannotated_calls: {type: 'File?', outputSource: manta/annotsv_unannotated_calls}
  mutect2_public_outputs: {type: 'File[]', outputSource: mutect2/mutect2_public_outputs}
  mutect2_protected_outputs: {type: 'File[]', outputSource: mutect2/mutect2_protected_outputs}
  mutect2_prepass_vcf: {type: 'File', outputSource: mutect2/mutect2_filtered_vcf}
  vardict_public_outputs: {type: 'File[]', outputSource: vardict/vardict_public_outputs}
  vardict_protected_outputs: {type: 'File[]', outputSource: vardict/vardict_protected_outputs}
  vardict_prepass_vcf: {type: 'File', outputSource: vardict/vardict_prepass_vcf}
  lancet_public_outputs: {type: 'File[]', outputSource: lancet/lancet_public_outputs}
  lancet_protected_outputs: {type: 'File[]', outputSource: lancet/lancet_protected_outputs}
  lancet_prepass_vcf: {type: 'File', outputSource: lancet/lancet_prepass_vcf}
  gatk_copy_ratio_segments_tumor: {type: 'File?', outputSource: gatk_cnv/called_copy_ratio_segments_tumor, doc: "Called copy-ratio-segments
      file. This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
      a row specifying the column headers contained in CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn, and the
      corresponding entry rows."}
  gatk_copy_ratio_segments_normal: {type: 'File?', outputSource: gatk_cnv/called_copy_ratio_segments_normal, doc: "Called copy-ratio-segments
      file. This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
      a row specifying the column headers contained in CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn, and the
      corresponding entry rows."}
  gatk_cnv_denoised_tumor_plot: {type: 'File?', outputSource: gatk_cnv/denoised_tumor_plot, doc: "Denoised-plot file that covers the
      entire range of the copy ratios"}
  gatk_cnv_denoised_normal_plot: {type: 'File?', outputSource: gatk_cnv/denoised_normal_plot, doc: "Denoised-plot file that covers
      the entire range of the copy ratios"}
  gatk_cnv_funcotated_called_file_tumor: {type: 'File?', outputSource: gatk_cnv/funcotated_called_file_tumor, doc: "TSV where each
      row is a segment and the annotations are the covered genes and which genes+exon is overlapped by the segment breakpoints."}
  gatk_cnv_funcotated_called_gene_list_file_tumor: {type: 'File?', outputSource: gatk_cnv/funcotated_called_gene_list_file_tumor,
    doc: "TSV where each row is a gene and the annotations are the covered genes and which genes+exon is overlapped by the segment
      breakpoints."}
steps:
  runtime_validator:
    run: ../tools/runtime_validator.cwl
    in:
      is_wgs:
        source: wgs_or_wxs
        valueFrom: |
          $(self == 'WGS')
      vardict: run_vardict
      mutect2: run_mutect2
      strelka2: run_strelka2
      preexisting_vcf:
        source: lancet_input_vcf
        valueFrom: $(self != null)
      lancet: run_lancet
      controlfreec: run_controlfreec
      cnvkit: run_cnvkit
      amplicon_architect: run_amplicon_architect
      theta2: run_theta2
      manta: run_manta
      gatk_cnv: run_gatk_cnv
      mosek_present:
        source: mosek_license_file
        valueFrom: |
          $(self != null)
      pon_present:
        source: count_panel_of_normals
        valueFrom: |
          $(self != null)
      exome_flag: exome_flag
      cnvkit_wgs_mode: cnvkit_wgs_mode
      i_flag: i_flag
      lancet_padding: lancet_padding
      lancet_window: lancet_window
      vardict_padding: vardict_padding
    out: [out_wgs, run_vardict, run_mutect2, run_strelka2, run_lancet, run_controlfreec, run_cnvkit, run_amplicon_architect, run_theta2,
      run_manta, run_gatk_cnv, out_exome_flag, out_cnvkit_wgs_mode, out_i_flag, out_lancet_padding, out_lancet_window, out_vardict_padding]
  samtools_cram2bam_plus_calmd_tumor:
    run: ../tools/samtools_calmd.cwl
    when: $(inputs.run_tool.some(function(e) { return e }) && (inputs.input_reads.basename.search(/.cram$/) != -1 || inputs.run_anyway) )
    in:
      run_tool:
        source: [runtime_validator/run_vardict, runtime_validator/run_lancet, runtime_validator/run_controlfreec, runtime_validator/run_cnvkit]
      run_anyway: run_calmd_bam
      input_reads: input_tumor_aligned
      threads:
        valueFrom: ${return 16;}
      reference: indexed_reference_fasta
    out: [bam_file]
  samtools_cram2bam_plus_calmd_normal:
    run: ../tools/samtools_calmd.cwl
    when: $(inputs.run_tool.some(function(e) { return e }) && (inputs.input_reads.basename.search(/.cram$/) != -1 || inputs.run_anyway) )
    in:
      run_tool:
        source: [runtime_validator/run_vardict, runtime_validator/run_lancet, runtime_validator/run_controlfreec, runtime_validator/run_cnvkit]
      run_anyway: run_calmd_bam
      input_reads: input_normal_aligned
      threads:
        valueFrom: ${return 16;}
      reference: indexed_reference_fasta
    out: [bam_file]
  prepare_regions_padded:
    run: ../sub_workflows/prepare_regions.cwl
    when: $(inputs.wgs_or_wxs == 'WXS' && inputs.run_tool.some(function(e) { return e }))
    in:
      wgs_or_wxs: wgs_or_wxs
      run_tool:
        source: [runtime_validator/run_mutect2, runtime_validator/run_strelka2, runtime_validator/run_vardict, runtime_validator/run_manta,
          runtime_validator/run_lancet]
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      calling_regions: calling_regions
      calling_padding:
        valueFrom: $(100)
      blacklist_regions: blacklist_regions
      break_bands_at_multiples_of:
        valueFrom: $(0)
      scatter_count:
        valueFrom: $(50)
    out: [prescatter_intervallist, prescatter_bed, prescatter_bedgz, scattered_intervallists, scattered_beds]
  prepare_regions_unpadded:
    run: ../sub_workflows/prepare_regions.cwl
    when: $(inputs.wgs_or_wxs == 'WGS' && inputs.run_tool.some(function(e) { return e }))
    in:
      wgs_or_wxs: wgs_or_wxs
      run_tool:
        source: [runtime_validator/run_mutect2, runtime_validator/run_strelka2, runtime_validator/run_vardict, runtime_validator/run_manta]
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      calling_regions: calling_regions
      blacklist_regions: blacklist_regions
      break_bands_at_multiples_of:
        valueFrom: $(80000000)
      scatter_count:
        valueFrom: $(50)
    out: [prescatter_intervallist, prescatter_bed, prescatter_bedgz, scattered_intervallists, scattered_beds]
  prepare_regions_unpadded_minibands:
    run: ../sub_workflows/prepare_regions.cwl
    when: $(inputs.wgs_or_wxs == 'WGS' && inputs.run_tool)
    in:
      wgs_or_wxs: wgs_or_wxs
      run_tool: runtime_validator/run_vardict
      calling_regions: prepare_regions_unpadded/prescatter_intervallist
      break_bands_at_multiples_of:
        valueFrom: $(20000)
      scatter_count:
        valueFrom: $(50)
    out: [prescatter_intervallist, prescatter_bed, prescatter_bedgz, scattered_intervallists, scattered_beds]
  prepare_regions_unpadded_cnv:
    run: ../sub_workflows/prepare_regions.cwl
    when: $(inputs.run_tool.some(function(e) { return e }))
    in:
      run_tool:
        source: [runtime_validator/run_controlfreec, runtime_validator/run_cnvkit]
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      calling_regions: calling_regions
      blacklist_regions: cnv_blacklist_regions
      scatter_count:
        valueFrom: $(0)
    out: [prescatter_intervallist, prescatter_bed, prescatter_bedgz, scattered_intervallists, scattered_beds]
  bedtools_intersect_germline:
    run: ../tools/bedtools_intersect.cwl
    when: $(inputs.run_tool.some(function(e) { return e }))
    in:
      run_tool:
        source: [runtime_validator/run_controlfreec, runtime_validator/run_cnvkit]
      input_vcf: b_allele
      output_basename: output_basename
      input_bed_file: prepare_regions_unpadded_cnv/prescatter_bed
      flag: runtime_validator/out_i_flag
    out: [intersected_vcf]
  gatk_filter_germline:
    run: ../tools/gatk_filter_germline_variant.cwl
    when: $(inputs.run_tool.some(function(e) { return e }))
    in:
      run_tool:
        source: [runtime_validator/run_controlfreec, runtime_validator/run_cnvkit]
      input_vcf: bedtools_intersect_germline/intersected_vcf
      reference_fasta: indexed_reference_fasta
      output_basename: output_basename
    out: [filtered_vcf, filtered_pass_vcf]
  vardict:
    run: ../sub_workflows/kfdrc_vardict_sub_wf.cwl
    when: $(inputs.run_vardict)
    in:
      run_vardict: runtime_validator/run_vardict
      indexed_reference_fasta: indexed_reference_fasta
      input_tumor_aligned:
        source: [ samtools_cram2bam_plus_calmd_tumor/bam_file, input_tumor_aligned]
        pickValue: first_non_null
      input_tumor_name: input_tumor_name
      old_tumor_name: old_tumor_name
      input_normal_aligned:
        source: [ samtools_cram2bam_plus_calmd_normal/bam_file, input_normal_aligned]
        pickValue: first_non_null
      input_normal_name: input_normal_name
      old_normal_name: old_normal_name
      output_basename: output_basename
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      bed_invtl_split:
        source: [prepare_regions_padded/scattered_beds, prepare_regions_unpadded_minibands/scattered_beds]
        valueFrom: |
          $(self.some(function(e) { return e != null }) ? self.filter(function(e) { return e != null })[0] : null)
      padding: runtime_validator/out_vardict_padding
      min_vaf: vardict_min_vaf
      select_vars_mode: select_vars_mode
      cpus: vardict_cpus
      ram: vardict_ram
      vep_cache: vep_cache
      vep_ram: vep_ram
      vep_cores: vep_cores
      vep_buffer_size: vep_buffer_size
      dbnsfp: dbnsfp
      dbnsfp_fields: dbnsfp_fields
      merged: merged
      cadd_indels: cadd_indels
      cadd_snvs: cadd_snvs
      run_cache_af: run_cache_af
      run_cache_existing: run_cache_existing
      retain_info: vardict_retain_info
      retain_fmt: vardict_retain_fmt
      retain_ann: vardict_retain_ann
      echtvar_anno_zips: echtvar_anno_zips
      bcftools_public_filter: bcftools_public_filter
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      disable_hotspot_annotation: disable_hotspot_annotation
      disable_vep_annotation: disable_vep_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      custom_enst: custom_enst
    out: [vardict_prepass_vcf, vardict_protected_outputs, vardict_public_outputs]
  mutect2:
    run: ../sub_workflows/kfdrc_mutect2_sub_wf.cwl
    when: $(inputs.run_mutect2)
    in:
      run_mutect2: runtime_validator/run_mutect2
      indexed_reference_fasta: indexed_reference_fasta
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      bed_invtl_split:
        source: [prepare_regions_padded/scattered_beds, prepare_regions_unpadded/scattered_beds]
        valueFrom: |
          $(self.some(function(e) { return e != null }) ? self.filter(function(e) { return e != null })[0] : null)
      af_only_gnomad_vcf: mutect2_af_only_gnomad_vcf
      exac_common_vcf: mutect2_exac_common_vcf
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      old_tumor_name: old_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      old_normal_name: old_normal_name
      exome_flag: runtime_validator/out_exome_flag
      output_basename: output_basename
      learnorientation_memory: learnorientation_memory
      getpileup_memory: getpileup_memory
      filtermutectcalls_memory: filtermutectcalls_memory
      select_vars_mode: select_vars_mode
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      vep_cache: vep_cache
      vep_ram: vep_ram
      vep_cores: vep_cores
      vep_buffer_size: vep_buffer_size
      dbnsfp: dbnsfp
      dbnsfp_fields: dbnsfp_fields
      merged: merged
      cadd_indels: cadd_indels
      cadd_snvs: cadd_snvs
      run_cache_af: run_cache_af
      run_cache_existing: run_cache_existing
      retain_info: mutect2_retain_info
      retain_fmt: mutect2_retain_fmt
      retain_ann: mutect2_retain_ann
      echtvar_anno_zips: echtvar_anno_zips
      bcftools_public_filter: bcftools_public_filter
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      custom_enst: custom_enst
      disable_vep_annotation: disable_vep_annotation
    out: [mutect2_filtered_stats, mutect2_filtered_vcf, mutect2_protected_outputs, mutect2_public_outputs]
  strelka2:
    run: ../sub_workflows/kfdrc_strelka2_sub_wf.cwl
    when: $(inputs.run_strelka2)
    in:
      run_strelka2: runtime_validator/run_strelka2
      indexed_reference_fasta: indexed_reference_fasta
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      hg38_strelka_bed:
        source: [prepare_regions_padded/prescatter_bedgz, prepare_regions_unpadded/prescatter_bedgz]
        valueFrom: |
          $(self.some(function(e) { return e != null }) ? self.filter(function(e) { return e != null })[0] : null)
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      manta_small_indels: manta/manta_small_indels
      use_manta_small_indels: use_manta_small_indels
      exome_flag: runtime_validator/out_exome_flag
      extra_arg: extra_arg
      strelka2_cores: strelka2_cores
      vep_cache: vep_cache
      vep_ram: vep_ram
      vep_cores: vep_cores
      vep_buffer_size: vep_buffer_size
      dbnsfp: dbnsfp
      dbnsfp_fields: dbnsfp_fields
      merged: merged
      cadd_indels: cadd_indels
      cadd_snvs: cadd_snvs
      run_cache_af: run_cache_af
      run_cache_existing: run_cache_existing
      retain_info: strelka2_retain_info
      retain_fmt: strelka2_retain_fmt
      retain_ann: strelka2_retain_ann
      output_basename: output_basename
      select_vars_mode: select_vars_mode
      echtvar_anno_zips: echtvar_anno_zips
      bcftools_public_filter: bcftools_public_filter
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      disable_hotspot_annotation: disable_hotspot_annotation
      disable_vep_annotation: disable_vep_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      custom_enst: custom_enst
    out: [strelka2_prepass_vcf, strelka2_protected_outputs, strelka2_public_outputs]
  prepare_regions_lancet_wgs:
    run: ../sub_workflows/prepare_regions.cwl
    when: $(inputs.wgs_or_wxs == 'WGS' && inputs.run_tool)
    in:
      wgs_or_wxs: wgs_or_wxs
      run_tool: runtime_validator/run_lancet
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      calling_regions: coding_sequence_regions
      supplement_vcfs:
        source: [strelka2/strelka2_protected_outputs, mutect2/mutect2_protected_outputs, lancet_input_vcf]
        valueFrom: |
          $(self.filter(function(e) { return e != null }).map(function(e) { return e.filter(function(i) { return i.basename.search(/(.vcf|.vcf.gz)$/) != -1 }) }).reduce(function(a,c) { return a.concat(c) }))
      blacklist_regions: blacklist_regions
      break_bands_at_multiples_of:
        valueFrom: $(0)
      scatter_count:
        valueFrom: $(50)
    out: [prescatter_intervallist, prescatter_bed, prescatter_bedgz, scattered_intervallists, scattered_beds]
  lancet:
    run: ../sub_workflows/kfdrc_lancet_sub_wf.cwl
    when: $(inputs.run_lancet)
    in:
      run_lancet: runtime_validator/run_lancet
      indexed_reference_fasta: indexed_reference_fasta
      input_tumor_aligned:
        source: [ samtools_cram2bam_plus_calmd_tumor/bam_file, input_tumor_aligned]
        pickValue: first_non_null
      input_tumor_name: input_tumor_name
      old_tumor_name: old_tumor_name
      input_normal_aligned:
        source: [ samtools_cram2bam_plus_calmd_normal/bam_file, input_normal_aligned]
        pickValue: first_non_null
      input_normal_name: input_normal_name
      old_normal_name: old_normal_name
      output_basename: output_basename
      select_vars_mode: select_vars_mode
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      bed_invtl_split:
        source: [prepare_regions_lancet_wgs/scattered_beds, prepare_regions_padded/scattered_beds]
        valueFrom: |
          $(self.some(function(e) { return e != null }) ? self.filter(function(e) { return e != null })[0] : null)
      ram: lancet_ram
      window: runtime_validator/out_lancet_window
      padding: runtime_validator/out_lancet_padding
      vep_cache: vep_cache
      vep_ram: vep_ram
      vep_cores: vep_cores
      vep_buffer_size: vep_buffer_size
      dbnsfp: dbnsfp
      dbnsfp_fields: dbnsfp_fields
      merged: merged
      cadd_indels: cadd_indels
      cadd_snvs: cadd_snvs
      run_cache_af: run_cache_af
      run_cache_existing: run_cache_existing
      retain_info: lancet_retain_info
      retain_fmt: lancet_retain_fmt
      retain_ann: lancet_retain_ann
      echtvar_anno_zips: echtvar_anno_zips
      bcftools_public_filter: bcftools_public_filter
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      disable_hotspot_annotation: disable_hotspot_annotation
      disable_vep_annotation: disable_vep_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      custom_enst: custom_enst
    out: [lancet_prepass_vcf, lancet_protected_outputs, lancet_public_outputs]
  controlfreec:
    run: ../sub_workflows/kfdrc_controlfreec_sub_wf.cwl
    when: $(inputs.run_controlfreec)
    in:
      run_controlfreec: runtime_validator/run_controlfreec
      input_tumor_aligned:
        source: [ samtools_cram2bam_plus_calmd_tumor/bam_file, input_tumor_aligned]
        pickValue: first_non_null
      input_tumor_name: input_tumor_name
      input_normal_aligned:
        source: [ samtools_cram2bam_plus_calmd_normal/bam_file, input_normal_aligned]
        pickValue: first_non_null
      wgs_or_wxs: wgs_or_wxs
      threads: cfree_threads
      output_basename: output_basename
      ploidy: cfree_ploidy
      mate_orientation_sample: cfree_mate_orientation_sample
      mate_orientation_control: cfree_mate_orientation_control
      calling_regions: prepare_regions_unpadded_cnv/prescatter_bed
      indexed_reference_fasta: indexed_reference_fasta
      b_allele: gatk_filter_germline/filtered_pass_vcf
      coeff_var: cfree_coeff_var
      contamination_adjustment: cfree_contamination_adjustment
      cfree_sex: cfree_sex
    out: [ctrlfreec_pval, ctrlfreec_config, ctrlfreec_pngs, ctrlfreec_bam_ratio, ctrlfreec_bam_seg, ctrlfreec_baf, ctrlfreec_info]
  cnvkit:
    run: ../sub_workflows/kfdrc_cnvkit_sub_wf.cwl
    when: $(inputs.run_cnvkit)
    in:
      run_cnvkit: runtime_validator/run_cnvkit
      input_tumor_aligned:
        source: [ samtools_cram2bam_plus_calmd_tumor/bam_file, input_tumor_aligned]
        pickValue: first_non_null
      tumor_sample_name: input_tumor_name
      input_normal_aligned:
        source: [ samtools_cram2bam_plus_calmd_normal/bam_file, input_normal_aligned]
        pickValue: first_non_null
      reference: indexed_reference_fasta
      normal_sample_name: input_normal_name
      capture_regions: prepare_regions_unpadded_cnv/prescatter_bed
      blacklist_regions: cnv_blacklist_regions
      wgs_mode: runtime_validator/out_cnvkit_wgs_mode
      b_allele_vcf: gatk_filter_germline/filtered_pass_vcf
      annotation_file: cnvkit_annotation_file
      output_basename: output_basename
      sex: cnvkit_sex
    out: [cnvkit_cnr, cnvkit_cnn_output, cnvkit_cns, cnvkit_calls, cnvkit_metrics, cnvkit_gainloss, cnvkit_seg, cnvkit_scatter_plot,
      cnvkit_diagram]
  amplicon_architect:
    run: ../workflow/kfdrc_production_amplicon_architect.cwl
    when: $(inputs.run_amplicon_architect)
    in:
      run_amplicon_architect: runtime_validator/run_amplicon_architect
      aa_data_repo: aa_data_repo
      aa_data_ref_version: aa_data_ref_version
      tumor_align_file: samtools_cram2bam_plus_calmd_tumor/bam_file
      output_basename: output_basename
      mosek_license_file: mosek_license_file
      reference: indexed_reference_fasta
      cnvkit_cns: cnvkit/cnvkit_cns
      male_input_flag:
        source: cnvkit_sex
        valueFrom: "$(self == 'y' ? true : null)"
      wgs_or_wxs: wgs_or_wxs
    out: [aa_cnv_seeds, aa_summary, aa_cycles, aa_graph, aa_sv_png, aa_classification_profiles, aa_gene_list]
  theta2_purity:
    run: ../sub_workflows/kfdrc_run_theta2_sub_wf.cwl
    when: $(inputs.run_theta2)
    in:
      run_theta2: runtime_validator/run_theta2
      tumor_cns: cnvkit/cnvkit_calls
      reference_cnn: cnvkit/cnvkit_cnn_output
      tumor_sample_name: input_tumor_name
      normal_sample_name: input_normal_name
      paired_vcf: vardict/vardict_prepass_vcf
      combined_include_expression: combined_include_expression
      combined_exclude_expression: combined_exclude_expression
      min_theta2_frac: min_theta2_frac
      output_basename: output_basename
    out: [theta2_adjusted_cns, theta2_adjusted_seg, theta2_subclonal_results, theta2_subclonal_cns, theta2_subclone_seg]
  expression_flatten_subclonal_results:
    run: ../tools/expression_flatten_file_list.cwl
    when: $(inputs.input_list != null)
    in:
      input_list: theta2_purity/theta2_subclonal_results
    out: [output]
  manta:
    run: ../sub_workflows/kfdrc_manta_sub_wf.cwl
    when: $(inputs.run_manta)
    in:
      run_manta: runtime_validator/run_manta
      indexed_reference_fasta: indexed_reference_fasta
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      hg38_strelka_bed:
        source: [prepare_regions_padded/prescatter_bedgz, prepare_regions_unpadded/prescatter_bedgz]
        valueFrom: |
          $(self.some(function(e) { return e != null }) ? self.filter(function(e) { return e != null })[0] : null)
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      old_tumor_name: old_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      old_normal_name: old_normal_name
      vep_cache: vep_cache
      output_basename: output_basename
      manta_memory: manta_memory
      manta_cores: manta_cores
      select_vars_mode: select_vars_mode
      annotsv_annotations_dir_tgz: annotsv_annotations_dir_tgz
    out: [manta_prepass_vcf, manta_pass_vcf, manta_small_indels, annotsv_annotated_calls, annotsv_unannotated_calls]
  gatk_cnv:
    run: ../sub_workflows/kfdrc_gatk_cnv_somatic_pair_wf.cwl
    when: $(inputs.run_gatk_cnv)
    in:
      run_gatk_cnv: runtime_validator/run_gatk_cnv
      input_aligned_reads_tumor: input_tumor_aligned
      input_aligned_reads_normal: input_normal_aligned
      reference_fasta: indexed_reference_fasta
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      input_interval_list: calling_regions
      input_exclude_interval_list: cnv_blacklist_regions
      bin_length:
        source: wgs_or_wxs
        valueFrom: |
          $(self == "WGS" ? 1000 : 0)
      common_sites: b_allele
      count_panel_of_normals: count_panel_of_normals
      output_basename: output_basename
      run_funcotatesegments: run_funcotatesegments
      funcotator_data_sources_tgz: funcotator_data_sources_tgz
      funcotator_minimum_segment_size: funcotator_minimum_segment_size
    out: [tumor_file_archive, modeled_segments_tumor, modeled_segments_tumor_plot, called_copy_ratio_segments_tumor, denoised_tumor_plot,
      normal_file_archive, modeled_segments_normal, modeled_segments_normal_plot, called_copy_ratio_segments_normal, denoised_normal_plot,
      funcotated_called_file_tumor, funcotated_called_gene_list_file_tumor]
$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 6
"sbg:license": Apache License 2.0
"sbg:publisher": KFDRC
"sbg:categories":
- AMPLICONARCHITECT
- BAM
- CNV
- CNVKIT
- CONTROLFREEC
- CRAM
- ECDNA
- FUNCOTATOR
- GATK
- INDEL
- LANCET
- MAF
- MANTA
- MUTECT2
- SNV
- SOMATIC
- STRELKA2
- SV
- THETA2
- VARDICT
- VCF
- VEP
"sbg:links":
- id: 'https://github.com/kids-first/kf-somatic-workflow/releases/tag/v5.1.0'
  label: github-release
