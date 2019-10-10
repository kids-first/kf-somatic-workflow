# KFDRC Somatic Whole Genome Sequence Analysis Workflow

![data service logo](https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcS9BnbvIsTkK3QlSGMDvlgu0tZQJ1q4crMvA-S3fcWfIq6y2d2Y)

This is the Kids First Data Resource Center (DRC) Whole Genome Sequencing (WGS) Somatic Workflow, which includes somatic variant calling, copy number variation (CNV), and structural variant (SV) calls. 
This workflow takes aligned cram input and performs somatic variant calling using Strelka2 and Mutect2, CNV estimation using ControlFreeC, and SV calls using Manta.
Somatic variant and SV call results are annoated using Variant Effect Predictor, with the Memorial Sloane Kettering Cancer Center (MSKCC) vcf2maf wrapper.

### Somatic Variant Calling:

[Strelka2](https://github.com/Illumina/strelka) v2.9.3 calls single nucelotide variants (SNVS) and insertions/deletions (INDELS).
[Mutect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php) v4.1.10 from the Broad institute also calls SNVS and INDELS.
Each caller has a different approach to variant calling, and together one can glean confident results.Strelka2 is run with default settings, similarly Mutect2 following Broad Best Practices, as of this [workflow](https://github.com/broadinstitute/gatk/blob/4.1.1.0/scripts/mutect2_wdl/mutect2.wdl).
Futhermore, each tool's results, in variant call format (vcf), are filtered on the `PASS` flag. 
The pre-`PASS` filtered results can still be obtained from the workflow in the eent the user wishes to keep some calls that failed `PASS` criteria.

### CNV Estimation:

[ControlFreeC](https://github.com/BoevaLab/FREEC) v11.6 is used for CNV estimation.
The tool portion of the workflow is a port from the Seven Bridges Genomics team, with a slight tweak in image outputs.
Also, the workflow wrapper limits what inputs and outputs are used based on our judgement of utility.
Outputs include raw ratio calls, copy number cals with p values assigned, b allele frequency data, as well as copy number and b allele frequency plots.

### SV Calls:

[Manta](https://github.com/Illumina/manta) v1.4.0 is used to call SVs. Output is also in vcf format, with calls filtered on `PASS`.
Default settings are used at run time.

### Variant Annotation

[Variant Effect Predictor](https://useast.ensembl.org/info/docs/tools/vep/index.html) release 93, wrapped by [vcf2maf](https://github.com/mskcc/vcf2maf) v1.6.16 is used to annotate somatic variant and SV calls.
Both the annotated vcf and maf file are made available.

### Tips To Run:

1) For input cram files, be sure to have indexed them beforehand as well.

2) For ControlFreeC, it is highly recommended that you supply a vcf file with germline calls, GATK Haplotype caller recommended.
Please also make sure the index for this file is available.
Also, a range of input ploidy possibilities for the inputs are needed.You can simply use `2`, or put in a range, as an array, like 2, 3, 4.

3) As a cavatica app, default references for hg38 are already pre-populated, as well as some default settings - i.e., number of threads, coefficient of variation input for ControlFreec, and `PASS` filter tool mode.

4) What is `select_vars_mode` you ask?On occasion, using GATK's `SelectVariants` tool will fail, so a simple `grep` mode on `PASS` can be used instead.
Related, `bcftools_filter_vcf` is built in as a convenience in case youe b allele frequency file has not been filtered on `PASS`.
You can use the `include_expression` `Filter="PASS"` to achieve this.

5) Suggested reference inputs are:

    - `reference_fasta`: [Homo_sapiens_assembly38.fasta](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK
    - `reference_dict`: [Homo_sapiens_assembly38.dict](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK
    - `wgs_calling_interval_list`: [wgs_calling_regions.hg38.interval_list](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK.*To create our 'wgs_canonical_calling_regions.hg38.interval_list', edit this file* by leaving only entries related to chr 1-22, X,Y, and M.M may need to be added.
    - `af_only_gnomad_vcf`: [af-only-gnomad.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/-gatk-best-practices/somatic-hg38) - need a valid google account, this is a link to the best practices google bucket from Broad GATK.
    - `exac_common_vcf`: [small_exac_common_3.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38) - need a valid google account, this is a link to the best practices google bucket from Broad GATK.
    - `hg38_strelka_bed`: [hg38_strelka.bed.gz'](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#extended-use-cases) - this link here has the bed-formatted text needed to copy to create this file.You will need to bgzip this file.
     - `threads`: 16
     - `chr_len`: hs38_chr.len, this a tsv file with chromosomes and their lengths.
      The first column must be chromosomes, optionally the secnod can be an alternate format of chromosomes.
      Last column must be chromosome length.
      Using the `hg38_strelka_bed`, and removing chrM can be a good source for this.
    - `coeff_var`: 0.05
    - `contamination_adjustment`: FALSE

6) Output files (Note, all vcf files that don't have an explicit index output have index files output as as secondary file.  In other words, they will be captured at the end of the workflow):

    - `strelka2_vep_vcf`: Variant effect predictor annotated vcf, filtered on `PASS`, somatic snv and indel call results from strelka2
    - `strelka2_vep_tbi`: Index file of above bgzipped vcf
    - `strelka2_prepass_vcf`: Somatic snv and indel call results with all `FILTER` categories for strelka2. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter.
    - `strelka2_vep_maf`: Mutation annotation file (maf) format of `strelka2_vep_vcf`
    - `mutect2_vep_vcf`: Variant effect predictor annotated vcf, filtered on `PASS`, somatic snv and indel call results from mutect2
    - `mutect2_vep_tbi`: Index file of above bgzipped vcf
    - `mutect2_prepass_vcf`: Somatic snv and indel call results with all `FILTER` categories for mutect2. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter.
    - `mutect2_vep_maf`: maf of format of `mutect2_vep_vcf`
    - `manta_vep_vcf`: Variant effect predictor annotated vcf, filtered on `PASS`, sv call results from manta
    - `manta_vep_tbi`: Index file of above bgzipped vcf
    - `manta_prepass_vcf`: SV results with all `FILTER` categories for manta. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter.
    - `manta_vep_maf`: maf of format of `manta_vep_vcf`
    - `ctrlfreec_pval`: CNV calls with copy number and p value confidence, a qualtitative "gain, loss, neutral" assignment, and genotype with uncertainty assigned from ControlFreeC.  See author manual for more details.
    - `ctrlfreec_config`: Config file used to run ControlFreeC.  Has some useful information on what parameters were used to run the tool.
    - `ctrlfreec_pngs`: Plots of b allele freqency (baf) log2 ratio and ratio of tumor/normal copy number coverage.  Pink line in the middle of ratio plots is the median ratio.
    - `ctrlfreec_bam_ratio`: Bam ratio text file.  Contain ratio, median ratio (used to inform `ctrlfreec_pval`), cnv estimates, baf estimate, and genotype estimate.
    - `ctrlfreec_bam_seg`: In-house generated seg file based on ratio file.  Provided asa a convenience for compatibility with tools that require inputs in legacy microarray format.
    - `ctrlfreec_baf`: baf estimations.
    - `ctrlfreec_info`: Contains useful run time information, like ploidy used for analysis, and window size

