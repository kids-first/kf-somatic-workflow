# KFDRC Somatic Whole Genome Sequence Analysis Workflow

![data service logo](https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcS9BnbvIsTkK3QlSGMDvlgu0tZQJ1q4crMvA-S3fcWfIq6y2d2Y)

This is the Kids First Data Resource Center (DRC) Whole Genome Sequencing (WGS) Somatic Workflow, which includes somatic variant calling, copy number variation (CNV), and structural variant (SV) calls. 
This workflow takes aligned cram input and performs somatic variant calling using Strelka2, Mutect2, Lancet, and VarDict Java, CNV estimation using ControlFreeC and CNVkit, and SV calls using Manta.
Somatic variant and SV call results are annotated using Variant Effect Predictor, with the Memorial Sloane Kettering Cancer Center (MSKCC) vcf2maf wrapper.
The `workflow/kfdrc_somatic_caller_workflow.cwl`, `workflow/kfdrc_cnvkit_plus_theta2.cwl`, `workflow/kfdrc_vardict_wf.cwl` and `workflow/kfdrc_lancet_wf.cwl` would run all tools described below for WGS.
We also have a whole [exome/targeted workflow](#KFDRC-Somatic-Whole Exome/Targeted-Sequence-Analysis-Workflow), which uses many of the same tools described for WGS.

### Somatic Variant Calling:

[Strelka2](https://github.com/Illumina/strelka) v2.9.3 calls single nucleotide variants (SNV) and insertions/deletions (INDEL).
[Mutect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php) v4.1.10 from the Broad institute calls SNV, multi-nucleotide variants (MNV, basically equal length substitutions with length > 1) and INDEL.
[Lancet](https://github.com/nygenome/lancet) v1.0.7 from the New York Genome Center (NYGC) calls SNV, MNV, and INDEL.
[VarDict Java](https://github.com/AstraZeneca-NGS/VarDictJava) v1.7.0 from AstraZeneca calls SNV, MNV, INDEL and more.
Each caller has a different approach to variant calling, and together one can glean confident results. Strelka2 is run with default settings, similarly Mutect2 following Broad Best Practices, as of this [workflow](https://github.com/broadinstitute/gatk/blob/4.1.1.0/scripts/mutect2_wdl/mutect2.wdl). Lancet is run in what I'd call an "Exome+" mode, based on the NYGC methods described [here](https://www.biorxiv.org/content/biorxiv/early/2019/04/30/623702.full.pdf). In short, regions from GENCODE gtf with feature annotations `exon`, `UTR`, and start/stop `codon` are used as intervals, as well as regions flanking hits from `strelka2` and `mutect2`. Lastly, VarDict Java run params follow the protocol that the [Blue Collar Bioinformatics](https://bcbio-nextgen.readthedocs.io/en/latest/index.html) uses, with the exception of using a min variant allele frequency (VAF) of 0.05 instead of 0.1, which we find to be relevant for rare cancer variant discovery. We also employ their false positive filtering methods.
Furthermore, each tool's results, in variant call format (vcf), are filtered on the `PASS` flag, with VarDict Java results additionally filtered for the flag `StrongSomatic`. Their results also include germline hits and other categories by default.
The pre-`PASS` filtered results can still be obtained from the workflow in the event the user wishes to keep some calls that failed `PASS` criteria.

### CNV Estimation:

[ControlFreeC](https://github.com/BoevaLab/FREEC) v11.6 is used for CNV estimation.
The tool portion of the workflow is a port from the [Seven Bridges Genomics](https://www.sevenbridges.com/) team, with a slight tweak in image outputs.
Also, the workflow wrapper limits what inputs and outputs are used based on our judgement of utility.
Outputs include raw ratio calls, copy number calls with p values assigned, b allele frequency data, as well as copy number and b allele frequency plots.
[CNVkit](https://cnvkit.readthedocs.io/en/stable/) v2.9.3 is a CNV second tool we currently use. [THeTa2](https://github.com/raphael-group/THetA) is used to inform and adjust copy number calls with purity estimations. For both tools, we take advantage of b allele frequency integration for copy number genotype estimation and increased CNV accuracy.

### SV Calls:

[Manta](https://github.com/Illumina/manta) v1.4.0 is used to call SVs. Output is also in vcf format, with calls filtered on `PASS`.
Default settings are used at run time.

### Variant Annotation

[Variant Effect Predictor](https://useast.ensembl.org/info/docs/tools/vep/index.html) release 93, wrapped by [vcf2maf](https://github.com/mskcc/vcf2maf) v1.6.17 is used to annotate somatic variant and SV calls.
Both the annotated vcf and maf file are made available.

### Tips To Run:

1) For input cram files, be sure to have indexed them beforehand as well.

2) For ControlFreeC, it is highly recommended that you supply a vcf file with germline calls, GATK Haplotype caller recommended.
Please also make sure the index for this file is available.
Also, a range of input ploidy possibilities for the inputs are needed. You can simply use `2`, or put in a range, as an array, like 2, 3, 4.
For mate orientation, you will need to specify, the drop down and tool doc explains your options.

3) As a cavatica app, default references for hg38 are already pre-populated, as well as some default settings - i.e., number of threads, coefficient of variation input for ControlFreec, and `PASS` filter tool mode.

4) What is `select_vars_mode` you ask? On occasion, using GATK's `SelectVariants` tool will fail, so a simple `grep` mode on `PASS` can be used instead.
Related, `bcftools_filter_vcf` is built in as a convenience in case your b allele frequency file has not been filtered on `PASS`.
You can use the `include_expression` `Filter="PASS"` to achieve this.

5) Suggested reference inputs are:

    - `reference_fasta`: [Homo_sapiens_assembly38.fasta](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK
    - `reference_dict`: [Homo_sapiens_assembly38.dict](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK
    - `annotation_file`: [refFlat_HG38.txt](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz) gunzip this filr from UCSC.  Needed for gene annotation in `CNVkit`
    - `wgs_calling_interval_list`: [wgs_calling_regions.hg38.interval_list](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK.*To create our 'wgs_canonical_calling_regions.hg38.interval_list', edit this file* by leaving only entries related to chr 1-22, X,Y, and M.M may need to be added.
    - `af_only_gnomad_vcf`: [af-only-gnomad.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/-gatk-best-practices/somatic-hg38) - need a valid google account, this is a link to the best practices google bucket from Broad GATK.
    - `exac_common_vcf`: [small_exac_common_3.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38) - need a valid google account, this is a link to the best practices google bucket from Broad GATK.
    - `hg38_strelka_bed`: [hg38_strelka.bed.gz'](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#extended-use-cases) - this link here has the bed-formatted text needed to copy to create this file. You will need to bgzip this file.
     - `vep_cache`: `homo_sapiens_vep_93_GRCh38.tar.gz` from ftp://ftp.ensembl.org/pub/release-93/variation/indexed_vep_cache/ - variant effect predictor cache.
     Current production workflow uses this version, and is compatible with the release used in the vcf2maf tool.
     - `threads`: 16
     - `chr_len`: hs38_chr.len, this a tsv file with chromosomes and their lengths. Should be limited to canonical chromosomes
      The first column must be chromosomes, optionally the second can be an alternate format of chromosomes.
      Last column must be chromosome length.
      Using the `hg38_strelka_bed`, and removing chrM can be a good source for this.
    - `coeff_var`: 0.05
    - `contamination_adjustment`: FALSE

6) Output files (Note, all vcf files that don't have an explicit index output have index files output as as secondary file.  In other words, they will be captured at the end of the workflow):

    - Simple variant callers
        - Strelka2:
            - `strelka2_vep_vcf`: Variant effect predictor annotated vcf, filtered on `PASS`, somatic snv and indel call results from strelka2
            - `strelka2_vep_tbi`: Index file of above bgzipped vcf
            - `strelka2_prepass_vcf`: Somatic snv and indel call results with all `FILTER` categories for strelka2. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter.
            - `strelka2_vep_maf`: Mutation annotation file (maf) format of `strelka2_vep_vcf`
        - Mutect2:
            - `mutect2_vep_vcf`: Variant effect predictor annotated vcf, filtered on `PASS`, somatic snv and indel call results from mutect2
            - `mutect2_vep_tbi`: Index file of above bgzipped vcf
            - `mutect2_prepass_vcf`: Somatic snv and indel call results with all `FILTER` categories for mutect2. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter.
            - `mutect2_vep_maf`: maf of format of `mutect2_vep_vcf`
        - VardictJava
            - `vardict_vep_somatic_only_vcf`: Variant effect predictor annotated vcf, filtered on `PASS` and `StrongSomatic` call results from VardictJava
            - `vardict_vep_somatic_only_tbi`: Index file of above bgzipped vcf
            - `vardict_vep_somatic_only_maf`: maf format of `vardict_vep_somatic_only_vcf`
            - `vardict_prepass_vcf`: All call results with all `FILTER` categories for VardictJava. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter and our `StrongSomatic` subset.
        - Lancet
          - `lancet_vep_vcf`: Variant effect predictor annotated vcf, filtered on `PASS`, somatic snv and indel call results from lancet
          - `lancet_vep_tbi`: Index file of above bgzipped vcf
          - `lancet_vep_maf`: maf format of `lancet_vep_vcf`
          - `lancet_prepass_vcf`: Somatic snv and indel call results with all `FILTER` categories for lancet. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter.
    - Structural variant callers
        - Manta
            - `manta_vep_vcf`: Variant effect predictor annotated vcf, filtered on `PASS`, sv call results from manta
            - `manta_vep_tbi`: Index file of above bgzipped vcf
            - `manta_prepass_vcf`: SV results with all `FILTER` categories for manta. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter.
            - `manta_vep_maf`: maf of format of `manta_vep_vcf`
    - Copy number variation callers
        - ControlFREEC
            - `ctrlfreec_pval`: CNV calls with copy number and p value confidence, a qualtitative "gain, loss, neutral" assignment, and genotype with uncertainty assigned from ControlFreeC.  See author manual for more details.
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


7) Docker images - the workflow tools will automatically pull them, but as a convenience are listed below:
    - `Strelka2`: obenauflab/strelka
    - `Mutect2` and all `GATK` tools: kfdrc/gatk:4.1.1.0
    - `Lancet`: kfdrc/lancet:1.0.7
    - `VarDict Java`: kfdrc/vardict:1.7.0
    - `ControlFreeC`: images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1
    - `CNVkit`: images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3
    - `THetA2`: kfdrc/theta2:0.7
    - `samtools`: kfdrc/samtools:1.9
    - `Variant Effect Predictor`: kfdrc/vep:r93_v2
    - `Manta`: kfdrc/manta:latest
    - `bcftools` and `vcftools`: kfdrc/bvcftools:latest 

# KFDRC Somatic Whole Exome/Targeted Sequence Analysis Workflow
Run this workflow if you have the exome/capture regions used.  If not, it is recommneded to use the WGS workflow instead, without the CNV and SV callers.  This is the Kids First Data Resource Center (DRC) Whole Exome/Targeted Sequencing (WXS) Somatic Workflow, which includes somatic variant calling and copy number variation (CNV). 
This workflow takes aligned cram input and performs somatic variant calling using Strelka2, Mutect2, Lancet, and VarDict Java, CNV estimation using ControlFreeC and CNVkit.
Somatic variant  call results are anntoated using Variant Effect Predictor, with the Memorial Sloane Kettering Cancer Center (MSKCC) vcf2maf wrapper.
The `workflow/kfdrc_production_WES_somatic_variant_cnv_wf.cwl` runs all tools described below.

### Somatic Variant Calling:

[Strelka2](https://github.com/Illumina/strelka) v2.9.3 calls single nucleotide variants (SNV) and insertions/deletions (INDEL).
[Mutect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php) v4.1.10 from the Broad institute calls SNV, multi-nucleotide variants (MNV, basically equal length substitutions with length > 1) and INDEL.
[Lancet](https://github.com/nygenome/lancet) v1.0.7 from the New York Genome Center (NYGC) calls SNV, MNV, and INDEL.
[VarDict Java](https://github.com/AstraZeneca-NGS/VarDictJava) v1.7.0 from AstraZeneca calls SNV, MNV, INDEL and more.
Each caller has a different approach to variant calling, and together one can glean confident results. Strelka2 is run with default settings, similarly Mutect2 following Broad Best Practices, as of this [workflow](https://github.com/broadinstitute/gatk/blob/4.1.1.0/scripts/mutect2_wdl/mutect2.wdl). Lancet is run with mostly default settings, with `max-indel-len=50`.  Padding and window size are adjustable, default 600 and 0 respectively, as input bed file is expected to be padded already. Lastly, VarDict Java run params follow the protocol that the [Blue Collar Informatics](https://bcbio-nextgen.readthedocs.io/en/latest/index.html) uses, with the exception of using a min variant allele frequency (VAF) of 0.05 instead of 0.1, which we find to be relevant for rare cancer variant discovery. We also employ their false positive filtering methods.
Furthermore, each tool's results, in variant call format (vcf), are filtered on the `PASS` flag, with VarDict Java results additionally filtered for the flag `StrongSomatic`. Their results also include germline hits and other categories by default.
The pre-`PASS` filtered results can still be obtained from the workflow in the event the user wishes to keep some calls that failed `PASS` criteria.

### CNV Estimation:

[ControlFreeC](https://github.com/BoevaLab/FREEC) v11.6 is used for CNV estimation.
The tool portion of the workflow is a port from the [Seven Bridges Genomics](https://www.sevenbridges.com/) team, with a slight tweak in image outputs.
Also, the workflow wrapper limits what inputs and outputs are used based on our judgement of utility.
Outputs include raw ratio calls, copy number calls with p values assigned, b allele frequency data, as well as copy number and b allele frequency plots.
[CNVkit](https://cnvkit.readthedocs.io/en/stable/) v2.9.3 is a CNV second tool we currently use. [THeTa2](https://github.com/raphael-group/THetA) is used to inform and adjust copy number calls with purity estimations. For both tools, we take advantage of b allele frequency integration for copy number genotype estimation and increased CNV accuracy.

### Variant Annotation

[Variant Effect Predictor](https://useast.ensembl.org/info/docs/tools/vep/index.html) release 93, wrapped by [vcf2maf](https://github.com/mskcc/vcf2maf) v1.6.17 is used to annotate somatic variant and SV calls.
Both the annotated vcf and maf file are made available.

### Tips To Run:

1) For input aligned files, cram or bam can be used. Be sure to index ahead of time.

2) For ControlFreeC, it is highly recommended that you supply a vcf file with germline calls, GATK Haplotype caller recommended.
Please also make sure the index for this file is available.
Also, a range of input ploidy possibilities for the inputs are needed. You can simply use `2`, or put in a range, as an array, like 2, 3, 4.
For mate orientation, you will need to specify, the drop down and tool doc explains your options.

3) As a cavatica app, default references for hg38 are already pre-populated, as well as some default settings - i.e., number of threads, coefficient of variation input for ControlFreec, and `PASS` filter tool mode.

4) What is `select_vars_mode` you ask? On occasion, using GATK's `SelectVariants` tool will fail, so a simple `grep` mode on `PASS` can be used instead.
Related, `bcftools_filter_vcf` is built in as a convenience in case your b allele frequency file has not been filtered on `PASS`.
You can use the `include_expression` `Filter="PASS"` to achieve this.

5) Suggested reference inputs are:

    - `reference_fasta`: [Homo_sapiens_assembly38.fasta](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK
    - `reference_dict`: [Homo_sapiens_assembly38.dict](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK
    - `annotation_file`: [refFlat_HG38.txt](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz) gunzip this filr from UCSC.  Needed for gene annotation in `CNVkit`
    - `padded_capture_regions`: Bed file with exome/targeted/capture regiond.  Recommend 100bp pad, for somatic variant calling
    - `unpadded_capture_regions`: Bed file with exome/targeted/capture regiond.  DO NOT pad, for copy number variation
    - `af_only_gnomad_vcf`: [af-only-gnomad.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/-gatk-best-practices/somatic-hg38) - need a valid google account, this is a link to the best practices google bucket from Broad GATK.
    - `exac_common_vcf`: [small_exac_common_3.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38) - need a valid google account, this is a link to the best practices google bucket from Broad GATK.
    - `hg38_strelka_bed`: This should be a bgzipped and tab-indexed version of `padded_capture_regions`
     - `vep_cache`: `homo_sapiens_vep_93_GRCh38.tar.gz` from ftp://ftp.ensembl.org/pub/release-93/variation/indexed_vep_cache/ - variant effect predictor cache.
     Current production workflow uses this version, and is compatible with the release used in the vcf2maf tool.
     - `threads`: 16
     - `cfree_chr_len`: hs38_chr.len, this a tsv file with chromosomes and their lengths. Should be limited to canonical chromosomes
      The first column must be chromosomes, optionally the second can be an alternate format of chromosomes.
      Last column must be chromosome length. This is needed for ControlFreeC
    - `coeff_var`: 0.05
    - `contamination_adjustment`: FALSE

6) Output files (Note, all vcf files that don't have an explicit index output have index files output as as secondary file.  In other words, they will be captured at the end of the workflow):

    - Simple variant callers
        - Strelka2:
            - `strelka2_vep_vcf`: Variant effect predictor annotated vcf, filtered on `PASS`, somatic snv and indel call results from strelka2
            - `strelka2_vep_tbi`: Index file of above bgzipped vcf
            - `strelka2_prepass_vcf`: Somatic snv and indel call results with all `FILTER` categories for strelka2. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter.
            - `strelka2_vep_maf`: Mutation annotation file (maf) format of `strelka2_vep_vcf`
        - Mutect2:
            - `mutect2_vep_vcf`: Variant effect predictor annotated vcf, filtered on `PASS`, somatic snv and indel call results from mutect2
            - `mutect2_vep_tbi`: Index file of above bgzipped vcf
            - `mutect2_prepass_vcf`: Somatic snv and indel call results with all `FILTER` categories for mutect2. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter.
            - `mutect2_vep_maf`: maf of format of `mutect2_vep_vcf`
        - VardictJava
            - `vardict_vep_somatic_only_vcf`: Variant effect predictor annotated vcf, filtered on `PASS` and `StrongSomatic` call results from VardictJava
            - `vardict_vep_somatic_only_tbi`: Index file of above bgzipped vcf
            - `vardict_vep_somatic_only_maf`: maf format of `vardict_vep_somatic_only_vcf`
            - `vardict_prepass_vcf`: All call results with all `FILTER` categories for VardictJava. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter and our `StrongSomatic` subset.
        - Lancet
          - `lancet_vep_vcf`: Variant effect predictor annotated vcf, filtered on `PASS`, somatic snv and indel call results from lancet
          - `lancet_vep_tbi`: Index file of above bgzipped vcf
          - `lancet_vep_maf`: maf format of `lancet_vep_vcf`
          - `lancet_prepass_vcf`: Somatic snv and indel call results with all `FILTER` categories for lancet. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter.
    - Copy number variation callers
        - ControlFreeeC
            - `ctrlfreec_pval`: CNV calls with copy number and p value confidence, a qualtitative "gain, loss, neutral" assignment, and genotype with uncertainty assigned from ControlFreeC.  See author manual for more details.
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

7) Docker images - the workflow tools will automatically pull them, but as a convenience are listed below:
    - `Strelka2`: obenauflab/strelka
    - `Mutect2` and all `GATK` tools: kfdrc/gatk:4.1.1.0
    - `Lancet`: kfdrc/lancet:1.0.7
    - `VarDict Java`: kfdrc/vardict:1.7.0
    - `ControlFreeC`: images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1
    - `CNVkit`: images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3
    - `samtools`: kfdrc/samtools:1.9
    - `Variant Effect Predictor`: kfdrc/vep:r93_v2
    - `bcftools` and `vcftools`: kfdrc/bvcftools:latest 
