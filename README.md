# Kids First DRC Somatic Variant Workflow

This repository contains the Kids First Data Resource Center (DRC) Somatic Variant Workflow, which includes somatic variant (SNV), copy number variation (CNV), and structural variant (SV) calls.
This workflow takes aligned cram input and performs somatic variant calling using Strelka2, Mutect2, Lancet, and VarDict Java, CNV estimation using ControlFreeC, CNVkit, and GATK, and SV calls using Manta.
For whole genome sequencing (WGS) data, the workflow will also predict extra chromosomal DNA (ecDNA) usiong AmpliconArchitect
Somatic variant call results are annotated with hotspots, assigned population frequencies using gnomAD AF, calculated gene models using Variant Effect Predictor, then added an additional MAF output using a modified version of MSKCCs vcf2maf.
See [annotation subworkflow doc](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc_annotation_subworkflow.md) for more details on annotation.

If you would like to run this workflow using the cavatica public app, a basic primer on running public apps can be found [here](https://www.notion.so/d3b/Starting-From-Scratch-Running-Cavatica-af5ebb78c38a4f3190e32e67b4ce12bb).
Alternatively, if you'd like to run it locally using `cwltool`, a basic primer on that can be found [here](https://www.notion.so/d3b/Starting-From-Scratch-Running-CWLtool-b8dbbde2dc7742e4aff290b0a878344d) and combined with app-specific info from the readme below.
This workflow is the current production workflow, equivalent to this [Cavatica public app](https://cavatica.sbgenomics.com/public/apps#cavatica/apps-publisher/kfdrc-somatic-variant-workflow)

![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

## Running WGS or WXS

The [combined workflow](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc-somatic-variant-workflow.cwl) is designed to be able to process either WGS or WXS inputs.
This functionality comes from usage of the `wgs_or_wxs` input enum. Depending on what is provided for this input, the tool will
set the appropriate default values and check that the user has provided the correct inputs. For example, if the user sets the
input to WGS the lancet_padding value will be defaulted to 300; alternatively, if the user sets the input to WXS the lancet_padding
value will be defaulted to 0. In either case, the user can override the defaults simply by providing their own value for lancet_padding
in the inputs.

The `wgs_or_wxs` flag also controls which inputs are used for certain steps. For example, the bed_interval input for Lancet comes
from different sources in the WGS and WXS pipelines. In the WGS pipeline separate processing is done ahead of time to generate
a new interval file. A tool in the workflow will take in the presumptive inputs for WGS and WXS modes. If the mode is WGS, then
the pipeline will pass on the file provided as the wgs_input and vice versa. If the wgs_input is missing and the mode is WGS, then
the pipeline will fail.

### WGS Run Fields

There are two WGS fields `wgs_calling_interval_list` and `lancet_calling_interval_bed`. If these are not provided in a WGS run,
the pipeline will fail.

### WXS Run Fields

There are two WXS fields `padded_capture_regions` and `unpadded_capture_regions`. If these are not provided in a WXS run,
the pipeline will fail.

### Standalone Somatic Workflows
Each tool used in the [combined workflow](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc-somatic-variant-workflow.cwl) can be run on its own. While the combined workflow calls all types of variant, each standalone caller only specializes in one class of variant.

| Workflow                                                                                                                                            | CNV | SNV | SV | ecDNA |
|-----------------------------------------------------------------------------------------------------------------------------------------------------|-----|-----|----|-------|
| [kfdrc-somatic-variant-workflow.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc-somatic-variant-workflow.cwl)     |  x  |  x  |  x |       |
| [kfdrc_production_cnvkit_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_cnvkit_wf.cwl)             |  x  |     |    |       |
| [kfdrc_production_controlfreec_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_controlfreec_wf.cwl) |  x  |     |    |       |
| [kfdrc_production_lancet_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_lancet_wf.cwl)             |     |  x  |    |       |
| [kfdrc_production_manta_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_manta_wf.cwl)               |     |     |  x |       |
| [kfdrc_production_mutect2_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_mutect2_wf.cwl)           |     |  x  |    |       |
| [kfdrc_production_strekla2_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_strekla2_wf.cwl)         |     |  x  |    |       |
| [kfdrc_production_theta2_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_theta2_wf.cwl)             |     |     |  x |       |
| [kfdrc_production_cnvkit_theta2_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_cnvkit_theta2_wf.cwl)|  x |     |    |       |
| [kfdrc_production_vardict_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_vardict_wf.cwl)           |     |  x  |    |       |
| [kfdrc_gatk_cnv_somatic_pair_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/sub_workflows/kfdrc_gatk_cnv_somatic_pair_wf.cwl)|  x  |     |    |       |
| [kfdrc_production_amplicon_architect.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_amplicon_architect.cwl)|     |     |      |   x  |
#### SNV Callers

- [Strelka2](https://github.com/Illumina/strelka/tree/v2.9.3) `v2.9.3`, from Illumina, calls single nucleotide variants (SNV) and insertions/deletions (INDEL)
  - See the [subworkflow doc](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc_strelka2_subworkflow.md) for more information
- [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360036730411-Mutect2) `v4.1.1.0`, from the Broad institute, calls SNV, multi-nucleotide variants (MNV, basically equal length substitutions with length > 1) and INDEL
  - This workflow will generate the interval lists needed to split up calling jobs to significantly reduce run time
  - Those intervals are used to run the [Mutect2 subworkflow](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc_mutect2_sub_wf.md)
- [Lancet](https://github.com/nygenome/lancet/releases/tag/v1.0.7) `v1.0.7`, from the New York Genome Center (NYGC), calls SNV, MNV, and INDEL
  - This workflow will generate the interval lists needed to split up calling jobs to significantly reduce run time
  - It will also convert cram input to bam input, if applicable
  - Intervals and bams are used as inputs to run the [Lancet subworkflow](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc_lancet_sub_wf.md)
- [VarDict Java](https://github.com/AstraZeneca-NGS/VarDictJava/tree/1.7.0) `v1.7.0`, from AstraZeneca, calls SNV, MNV, INDEL and more
  - This workflow will generate the interval lists needed to split up calling jobs to significantly reduce run time
  - It will also convert cram input to bam input, if applicable
  - Intervals and bams are used as inputs to run the [VarDict Java subworkflow](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc_vardict_sub_wf.md)

Each caller has a different approach to variant calling, and together one can glean confident results.
**After running this overall workflow, we recommend running our [consensus calling workflow](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc-consensus-calling.md) for a balance of sensitivity and specificity overall.**

#### CNV Estimators

- [ControlFreeC](https://github.com/BoevaLab/FREEC) `v11.6` is used for CNV estimation.
The tool portion of the workflow is a port from the [Seven Bridges Genomics](https://www.sevenbridges.com/) team, with a slight tweak in image outputs.
Also, the workflow wrapper limits what inputs and outputs are used based on our judgement of utility.
Outputs include raw ratio calls, copy number calls with p values assigned, b allele frequency data, as well as copy number and b allele frequency plots.
- [CNVkit](https://cnvkit.readthedocs.io/en/v0.9.3/) `v0.9.3` is a CNV second tool we currently use.
- [THeTa2](https://github.com/kids-first/THetA/tree/v0.7.1) is used to inform and adjust copy number calls from CNVkit with purity estimations.
- [GATK CNV](https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants) uses GATK 4.2.4.1 to call somatic CNVs using a Panel of Normals created using [this workflow](https://github.com/kids-first/kf-gatk-cnv-wf/blob/master/workflows/kf_create_cnv_pon_wf.cwl). **Note: If a PON is not provided, GATK CNV will be skipped!**

For ControlFreeC and CNVkit, we take advantage of b allele frequency (using the gVCF created by our [alignment and haplotypecaller workflows](https://github.com/kids-first/kf-alignment-workflow), **then** running our [Single Sample Genotyping Workflow](https://cavatica.sbgenomics.com/public/apps#cavatica/apps-publisher/kfdrc-single-sample-genotyping-wf/) on the gVCF) integration for copy number genotype estimation and increased CNV accuracy. Additionally these tools make use of the `unpadded_capture_regions` to provide the canonical calling regions.

#### ecDNA Prediction
 - [AmpliconArchitect](https://github.com/jluebeck/PrepareAA/blob/master/GUIDE.md) is used to predict ecDNAs
 - Workflfow documentation available [here](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc_amplicon_architect.md)

#### SV Callers and Annotators

- [Manta](https://github.com/Illumina/manta/tree/v1.4.0) `v1.4.0` is used to call SVs. Output is also in vcf format, with calls filtered on `PASS`.
Default settings are used at run time.
- [AnnotSV](https://github.com/lgmgeo/AnnotSV/releases/tag/v3.1.1) `v3.1.1` is used to annotate the calls in the Manta Prepass VCF.

#### Variant Annotation
Please see the [annotation workflow doc](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc_annotation_wf.md).
Both the annotated vcf and maf file are made available.

### Tips to Run:

1. For input cram files, be sure to have indexed them beforehand

1. When in doubt, all of our reference files can be obtained from here: https://cavatica.sbgenomics.com/u/kfdrc-harmonization/kf-references/

1. For ControlFreeC, it is highly recommended that you supply a vcf file with germline calls, GATK Haplotype caller recommended.
Please also make sure the index for this file is available.
Also, a range of input ploidy possibilities for the inputs are needed. You can simply use `2`, or put in a range, as an array, like 2, 3, 4.
For mate orientation, you will need to specify, the drop down and tool doc explains your options.

1. As a cavatica app, default references for hg38 are already pre-populated, as well as some default settings - i.e., number of threads, coefficient of variation input for ControlFreec, and `PASS` filter tool mode.

1. What is `select_vars_mode` you ask? On occasion, using GATK's `SelectVariants` tool will fail, so a simple `grep` mode on `PASS` can be used instead.
Related, `bcftools_filter_vcf` is built in as a convenience in case your b allele frequency file has not been filtered on `PASS`.
You can use the `include_expression` `Filter="PASS"` to achieve this.

1. The `SM:` sample IDs in the input alignment files **must** match `input_tumor_name` and `input_normal_name`. If you need the output VCF sample names to differ from the input, enter what is currently in `SM:` into the `old_tumor_name` and `old_normal_name` fields, and enter the desired sample names in `input_tumor_name` and `input_normal_name`. The workflow will rename the samples in the VCF for you.

1. Again, when in doubt our reference inputs can be obtained from [here](https://cavatica.sbgenomics.com/u/kfdrc-harmonization/kf-references/). Suggested reference inputs are:

    - `reference_fasta`: [Homo_sapiens_assembly38.fasta](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK
    - `reference_dict`: [Homo_sapiens_assembly38.dict](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK
    - `annotation_file`: [refFlat_HG38.txt](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz) gunzip this file from UCSC.  Needed for gene annotation in `CNVkit`
    - `wgs_calling_interval_list`: [wgs_calling_regions.hg38.interval_list](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK.*To create our 'wgs_canonical_calling_regions.hg38.interval_list', edit this file* by leaving only entries related to chr 1-22, X,Y, and M.M may need to be added.
    - `lancet_calling_interval_bed`: `GRCh38.gencode.v31.CDS.merged.bed`.  As described at the beginning, for WGS, it's highly recommended to use CDS bed, and supplement with region calls from Strelka2 & Mutect2. Our reference was obtained from GENCODE, [release 31](https://www.gencodegenes.org/human/release_31.html) using this gtf file [gencode.v31.primary_assembly.annotation.gtf.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.primary_assembly.annotation.gtf.gz) and parsing features for `UTR`, `start codon`, `stop codon`, and `exon`, then using bedtools sort and merge after converting coordinates into bed format.
    - `af_only_gnomad_vcf`: [af-only-gnomad.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/-gatk-best-practices/somatic-hg38) - need a valid google account, this is a link to the best practices google bucket from Broad GATK.
    - `exac_common_vcf`: [small_exac_common_3.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38) - need a valid google account, this is a link to the best practices google bucket from Broad GATK.
    - `hg38_strelka_bed`: [hg38_strelka.bed.gz'](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#extended-use-cases) - this link here has the bed-formatted text needed to copy to create this file. You will need to bgzip this file.
    - `extra_arg`: This can be used to add special params to strelka2. It is currently more of an "unsticking param". For edge cases where strelka2 seems to hang, setting this to `--max-input-depth 10000` can balance performance and consistency in results
    - strelka2_cores: `18`. This default is already set, but can be changed if desired.
    - `vep_cache`: `homo_sapiens_vep_93_GRCh38.tar.gz` from ftp://ftp.ensembl.org/pub/release-93/variation/indexed_vep_cache/ - variant effect predictor cache.
     Current production workflow uses this version.
    - `threads`: 16
    - `chr_len`: hs38_chr.len, this a tsv file with chromosomes and their lengths. Should be limited to canonical chromosomes
      The first column must be chromosomes, optionally the second can be an alternate format of chromosomes.
      Last column must be chromosome length.
      Using the `hg38_strelka_bed`, and removing chrM can be a good source for this.
    - `coeff_var`: 0.05
    - `contamination_adjustment`: FALSE
    - `genomic_hotspots`: `tert.bed`. Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots. This can be obtained from our cavatica reference project
    - `protein_snv_hotspots`: [hotspots_v2.xls](https://www.cancerhotspots.org/files/hotspots_v2.xls). Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots. Recommend pulling the two relevant columns for SNVs only, and convert to tsv
    - `protein_indel_hotspots`: [hotspots_v2.xls](https://www.cancerhotspots.org/files/hotspots_v2.xls). Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspotsRecommend pulling the two relevant columns for SNVs only, and convert to tsv
    bcftools_annot_columns: `"INFO/gnomad_3_1_1_AC:=INFO/AC,INFO/gnomad_3_1_1_AN:=INFO/AN,INFO/gnomad_3_1_1_AF:=INFO/AF,INFO/gnomad_3_1_1_nhomalt:=INFO/nhomalt,INFO/gnomad_3_1_1_AC_popmax:=INFO/AC_popmax,INFO/gnomad_3_1_1_AN_popmax:=INFO/AN_popmax,INFO/gnomad_3_1_1_AF_popmax:=INFO/AF_popmax,INFO/gnomad_3_1_1_nhomalt_popmax:=INFO/nhomalt_popmax,INFO/gnomad_3_1_1_AC_controls_and_biobanks:=INFO/AC_controls_and_biobanks,INFO/gnomad_3_1_1_AN_controls_and_biobanks:=INFO/AN_controls_and_biobanks,INFO/gnomad_3_1_1_AF_controls_and_biobanks:=INFO/AF_controls_and_biobanks,INFO/gnomad_3_1_1_AF_non_cancer:=INFO/AF_non_cancer,INFO/gnomad_3_1_1_primate_ai_score:=INFO/primate_ai_score,INFO/gnomad_3_1_1_splice_ai_consequence:=INFO/splice_ai_consequence"`. csv string of columns from annotation to port into the input vcf
    - `bcftools_annot_vcf`: `gnomad_3.1.1.vwb_subset.vcf.gz`. An export of the gnomAD v3.1.1 genomes reference made available from the KD workbench.
    - `bcftools_public_filter`: `'FILTER="PASS"|INFO/HotSpotAllele=1'`. This phrase will allow `PASS` only **or** `HotSpotAllele` variants into the public version of variant call output.
    - `gatk_filter_name`: `["NORM_DP_LOW", "GNOMAD_AF_HIGH"]`. These correspond to the recommended filter expression
    - `gatk_filter_expression`: `["vc.getGenotype('<input_normal_name> ').getDP() <= 7"), "gnomad_3_1_1_AF > 0.001"]`. Array of filter expressions to establish criteria to tag variants with. See [annotation subworkflow docs](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc_annotation_subworkflow.md) for a more detailed explanation. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration for general JEXL syntax
    - `disable_hotspot_annotation`: false
    - `maf_center`: `"."`. Sequencing center of variant called
    - `funcotator_data_sources_tgz`: [funcotator_dataSources.v1.6.20190124s.tar.gz](https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator) - need a valid google account, this is a link to pre-packaged datasources from the Broad Institute. Any of the tar.gz files will do.
    - `annotsv_annotations_dir_tgz`: [annotsv_311_annotations_dir.tgz] - These annotations are simply those from the install-human-annotation installation process run during AnnotSV installation. Specifically these are the annotations installed with v3.1.1 of the software. Newer or older annotations can be slotted in here as needed.
    - `aa_data_repo`: [GRCh38_indexed.tar.gz](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/) - amplicon architect reference tar ball. Only needed/used if input is WGS
    - `mosek_license_file`: [mosek.lic](https://www.mosek.com/license/request/) - required if amplicon architect is run. License is good for one year and is renewable
    - `aa_data_ref_version`: `"GRCh38"`. Genome reference version used

1. Output files (Note, all vcf files that don't have an explicit index output have index files output as as secondary file.  In other words, they will be captured at the end of the workflow):

    - Simple variant callers
        - Strelka2:
            - `strelka2_prepass_vcf`: Combined SNV + INDEL file with renamed Sample IDs. Has all soft `FILTER` values generated by variant caller
            - `strelka2_protected_outputs`: Array of files containing MAF format of PASS hits, `PASS` VCF with annotation pipeline soft `FILTER`-added values, and VCF index
            - `strelka2_public_outputs`: Same as above, except MAF and VCF have had entries with soft `FILTER` values removed
        - Mutect2:
            - `mutect2_filtered_vcf`: VCF with SNV, MNV, and INDEL variant calls. Contains all soft `FILTER` values generated by variant caller. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter
            - `mutect2_protected_outputs`: Array of files containing MAF format of PASS hits, `PASS` VCF with annotation pipeline soft `FILTER`-added values, and VCF index
            - `mutect2_public_outputs`: Same as above, except MAF and VCF have had entries with soft `FILTER` values removed
        - VardictJava
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
            - `manta_pass_tbi`: Index file of above bgzipped vcf
            - `manta_prepass_vcf`: SV results with all `FILTER` categories for manta. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter
        - AnnotSV
            - `annotsv_annotated_calls`: This file contains all records from the `manta_pass_vcf` that AnnotSV could annotate.
            - `annotsv_unannotated_calls`: This file contains all records from the `manta_pass_vcf` that AnnotSV could not annotate.
    - Copy number variation callers
        - ControlFREEC
            - `ctrlfreec_pval`: CNV calls with copy number and p value confidence, a qualitative "gain, loss, neutral" assignment, and genotype with uncertainty assigned from ControlFreeC.  See author manual for more details.
            - `ctrlfreec_config`: Config file used to run ControlFreeC.  Has some useful information on what parameters were used to run the tool.
            - `ctrlfreec_pngs`: Plots of b allele frequency (baf) log2 ratio and ratio of tumor/normal copy number coverage.  Pink line in the middle of ratio plots is the median ratio.
            - `ctrlfreec_bam_ratio`: Bam ratio text file.  Contain ratio, median ratio (used to inform `ctrlfreec_pval`), cnv estimates, baf estimate, and genotype estimate.
            - `ctrlfreec_bam_seg`: In-house generated seg file based on ratio file.  Provided asa a convenience for compatibility with tools that require inputs in legacy microarray format.
            - `ctrlfreec_baf`: baf estimations.
            - `ctrlfreec_info`: Contains useful run time information, like ploidy used for analysis, and window size
        - CNVkit
            - `cnvkit_cnr`: Copy number ratio
            - `cnvkit_cnn_output`: Normal/control sample copy number
            - `cnvkit_calls`: Tumor/sample copy number
            - `cnvkit_metrics`: Basic seg count and call stats
            - `cnvkit_gainloss`: Per-gene log2 ratio
            - `cnvkit_seg`: Classic microarray-style seg file
            - `theta2_calls`: Purity-adjusted CNVkit copy number calls based on theta2 results
            - `theta2_seg`: Purity-adjusted CNVkit seg file based on theta results
            - `theta2_subclonal_results`: Theta2 Subclone purity results
            - `theta2_subclonal_cns`: Theta2 sublone cns
            - `theta2_subclone_seg`: Theta subclone seg file
        - GATK CNV
            - `gatk_copy_ratio_segments_tumor`: Called copy-ratio-segments file for the tumor sample. This is a tab-separated values (TSV) file with SAM-style header containing a read group sample name, a sequence dictionary, a row specifying the column headers output by `CallCopyRatioSegments` and the corresponding entry rows.
            - `gatk_copy_ratio_segments_normal`: Called copy-ratio-segments file for the normal sample. This is a tab-separated values (TSV) file with SAM-style header containing a read group sample name, a sequence dictionary, a row specifying the column headers output by `CallCopyRatioSegments` and the corresponding entry rows.
            - `gatk_cnv_denoised_tumor_plot`: Denoised-plot file that covers the entire range of the copy ratios
            - `gatk_cnv_denoised_normal_plot`: Denoised-plot file that covers the entire range of the copy ratios
            - `gatk_cnv_funcotated_called_file_tumor`: TSV where each row is a segment and the annotations are the covered genes and which genes+exon is overlapped by the segment breakpoints.
            - `gatk_cnv_funcotated_called_gene_list_file_tumor`: TSV where each row is a gene and the annotations are the covered genes and which genes+exon is overlapped by the segment breakpoints.
    - ecDNA Prediction
        - `aa_summary`: Summary for all amplicons detected by AA
        - `aa_cycles`: Text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant source) and their copy counts
        - `aa_graph`: A text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts
        - `aa_sv_png`: PNG image file displaying the SV view of AA
        - `aa_classification_profiles`: abstract classification of the amplicon

1. Docker images - the workflow tools will automatically pull them, but as a convenience are listed below:
    - `Strelka2`: pgc-images.sbgenomics.com/d3b-bixu/strelka
    - `add common fields to Strelka2`: pgc-images.sbgenomics.com/d3b-bixu/add-strelka2-fields:1.0.0
    - `Mutect2` and `GATK SNV` tools: pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
    - `GATK CNV`: broadinstitute/gatk:4.2.4.1
    - `Lancet`: pgc-images.sbgenomics.com/d3b-bixu/lancet:1.0.7
    - `VarDict Java`: pgc-images.sbgenomics.com/d3b-bixu/vardict:1.7.0
    - `ControlFreeC`: images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1
    - `CNVkit`: images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3
    - `THetA2`: pgc-images.sbgenomics.com/d3b-bixu/theta2:0.7.1
    - `samtools`: pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9
    - `Variant Effect Predictor`: pgc-images.sbgenomics.com/d3b-bixu/vep:r93.7
    - `Kids Fist VCF2MAF`: pgc-images.sbgenomics.com/d3b-bixu/kf_vcf2maf:v1.0.3
    - `Manta`: pgc-images.sbgenomics.com/d3b-bixu/manta:1.4.0
    - `bcftools` and `vcftools`: pgc-images.sbgenomics.com/d3b-bixu/bvcftools:latest
    - `annotsv`: pgc-images.sbgenomics.com/d3b-bixu/annotsv:3.1.1
    - `AmpliconArchitect`: jluebeck/prepareaa:v0.1203.10

1. For highly complex samples, some tools have shown themselves to require memory allocation adjustments:
   Manta, GATK LearnReadOrientationModel, GATK GetPileupSummaries, GATK FilterMutectCalls and Vardict.
   Optional inputs exist to expand the memory allocation for these jobs: manta_memory, learnorientation_memory,
   getpileup_memory, filtermutectcalls_memory, and vardict_ram, respectively.
   For the java tools (Vardict and GATK), these values represent limits on the memory that can
   be used for the respective jobs. Tasks will go to these values and not exceed it. They may succeed or fail,
   but they will not exceed the limit established. The memory allocations for these is hardcapped. The memory
   allocation option for Manta, conversely, is a soft cap. The memory requested will be allocated for the job
   on a particular machine but once the task is running on the machine it may exceed that requested value. For example,
   if Manta's memory allocation is set to 10 GB it will have 10 GB allocated to it at task creation, but, if the
   task ends up running on a machine with more memory available, the task may use it. Setting a value here for Manta
   will not prevent Manta from taking more than that value. The memory usage in Manta is limited by the machine hardware.
   As such the option for Manta memory allocation is described as soft cap. For more information on Manta resource
   usage see their [documentation](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#runtime-hardware-requirements).

1. The optional `b_allele` file can be generated using our [Single Sample Genotyping Workflow](https://cavatica.sbgenomics.com/public/apps#cavatica/apps-publisher/kfdrc-single-sample-genotyping-wf/).

## Other Resources
- dockerfiles: https://github.com/d3b-center/bixtools