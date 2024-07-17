# Kids First DRC Somatic Variant Workflow

<p align="center">
  <img src="https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png">
</p>

This repository contains the Kids First Data Resource Center (DRC) Somatic Variant Workflow, which includes somatic variant (SNV), copy number variation (CNV), and structural variant (SV) calls.
Benchmarking of our SNV callers and [consensus methods](docs/kfdrc-consensus-calling.md) can be found [here](docs/SOMATIC_SNV_BENCHMARK.md).
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
   - `protein_snv_hotspots`: `kfdrc_protein_snv_cancer_hotspots_20240718.txt`. Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots. File header contains generation history
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
   - `gatk_filter_expression`: `["vc.getGenotype('<input_normal_name> ').getDP() <= 7"), "gnomad_3_1_1_AF != '.' && gnomad_3_1_1_AF > 0.001 && gnomad_3_1_1_FILTER=='PASS'"]`. Array of filter expressions to establish criteria to tag variants with. See [annotation subworkflow docs](./docs/kfdrc_annotation_subworkflow.md) for a more detailed explanation. For more information on filter expressions, see the [GATK JEXL docs](https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration).
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
