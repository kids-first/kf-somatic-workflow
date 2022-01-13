cwlVersion: v1.0
class: Workflow
id: kfdrc-somatic-variant-workflow
label: Kids First DRC Somatic Variant Workflow
doc: |
  # Kids First DRC Somatic Variant Workflow

  This repository contains the Kids First Data Resource Center (DRC) Somatic Variant Workflow, which includes somatic variant (SNV), copy number variation (CNV), and structural variant (SV) calls.
  This workflow takes aligned cram input and performs somatic variant calling using Strelka2, Mutect2, Lancet, and VarDict Java, CNV estimation using ControlFreeC and CNVkit, and SV calls using Manta.
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

  ### WGS-only Fields

  There are two WGS only fields `wgs_calling_interval_list` and `lancet_calling_interval_bed`. If these are not provided in a WGS run,
  the pipeline will fail.

  ### WXS-only Fields

  There are two WXS only fields `padded_capture_regions` and `unpadded_capture_regions`. If these are not provided in a WXS run,
  the pipeline will fail.

  ### Standalone Somatic Workflows
  Each tool used in the [combined workflow](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc-somatic-variant-workflow.cwl) can be run on its own. While the combined workflow calls all types of variant, each standalone caller only specializes in one class of variant.

  | Workflow                                                                                                                                            | CNV | SNV | SV |
  |-----------------------------------------------------------------------------------------------------------------------------------------------------|-----|-----|----|
  | [kfdrc-somatic-variant-workflow.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc-somatic-variant-workflow.cwl)     |  x  |  x  |  x |
  | [kfdrc_production_cnvkit_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_cnvkit_wf.cwl)             |  x  |     |    |
  | [kfdrc_production_controlfreec_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_controlfreec_wf.cwl) |  x  |     |    |
  | [kfdrc_production_lancet_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_lancet_wf.cwl)             |     |  x  |    |
  | [kfdrc_production_manta_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_manta_wf.cwl)               |     |     |  x |
  | [kfdrc_production_mutect2_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_mutect2_wf.cwl)           |     |  x  |    |
  | [kfdrc_production_strekla2_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_strekla2_wf.cwl)         |     |  x  |    |
  | [kfdrc_production_theta2_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_theta2_wf.cwl)             |     |     |  x |
  | [kfdrc_production_cnvkit_theta2_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_cnvkit_theta2_wf.cwl)             |   x |     |   |
  | [kfdrc_production_vardict_wf.cwl](https://github.com/kids-first/kf-somatic-workflow/blob/master/workflow/kfdrc_production_vardict_wf.cwl)           |     |  x  |    |

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

  For ControlFreeC and CNVkit, we take advantage of b allele frequency (from the gVCF created by our [alignment and haplotypecaller workflows](https://github.com/kids-first/kf-alignment-workflow)) integration for copy number genotype estimation and increased CNV accuracy.

  #### SV Callers

  - [Manta](https://github.com/Illumina/manta/tree/v1.4.0) `v1.4.0` is used to call SVs. Output is also in vcf format, with calls filtered on `PASS`.
  Default settings are used at run time.

  #### Variant Annotation
  Please see the [annotation subworkflow doc](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc_annotation_subworkflow.md).
  Both the annotated vcf and maf file are made available.

  ### Tips to Run:

  1. For input cram files, be sure to have indexed them beforehand as well.

  1. When in doubt, all of our reference files can be obtained from here: https://cavatica.sbgenomics.com/u/kfdrc-harmonization/kf-references/

  1. For ControlFreeC, it is highly recommended that you supply a vcf file with germline calls, GATK Haplotype caller recommended.
  Please also make sure the index for this file is available.
  Also, a range of input ploidy possibilities for the inputs are needed. You can simply use `2`, or put in a range, as an array, like 2, 3, 4.
  For mate orientation, you will need to specify, the drop down and tool doc explains your options.

  1. As a cavatica app, default references for hg38 are already pre-populated, as well as some default settings - i.e., number of threads, coefficient of variation input for ControlFreec, and `PASS` filter tool mode.

  1. What is `select_vars_mode` you ask? On occasion, using GATK's `SelectVariants` tool will fail, so a simple `grep` mode on `PASS` can be used instead.
  Related, `bcftools_filter_vcf` is built in as a convenience in case your b allele frequency file has not been filtered on `PASS`.
  You can use the `include_expression` `Filter="PASS"` to achieve this.

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
      `protein_indel_hotspots`: [hotspots_v2.xls](https://www.cancerhotspots.org/files/hotspots_v2.xls). Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspotsRecommend pulling the two relevant columns for SNVs only, and convert to tsv
      bcftools_annot_columns: `"INFO/AF"`. csv string of columns from annotation to port into the input vcf
      `bcftools_annot_vcf`: `af-only-gnomad.hg38.vcf.gz`. Yes, same as the Mutect2 reference listed above.
      `bcftools_public_filter`: `'FILTER="PASS"|INFO/HotSpotAllele=1'`. This phrase will allow `PASS` only **or** `HotSpotAllele` variants into the public version of variant call output.
      `gatk_filter_name`: `["NORM_DP_LOW", "GNOMAD_AF_HIGH"]`. These correspond to the recommended filter expression
      `gatk_filter_expression`: `["vc.getGenotype('<input_normal_name> ').getDP() <= 7"), "AF > 0.001"]`. Array of filter expressions to establish criteria to tag variants with. See [annotation subworkflow docs](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc_annotation_subworkflow.md) for a more detailed explanation. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration for general JEXL syntax
      `disable_hotspot_annotation`: false
      `maf_center`: `"."`. Sequencing center of variant called


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
              - `manta_vep_vcf`: SV call filtered on `PASS`, from manta
              - `manta_vep_tbi`: Index file of above bgzipped vcf
              - `manta_prepass_vcf`: SV results with all `FILTER` categories for manta. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter
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

  1. Docker images - the workflow tools will automatically pull them, but as a convenience are listed below:
      - `Strelka2`: pgc-images.sbgenomics.com/d3b-bixu/strelka
      - `add common fields to Strelka2`: pgc-images.sbgenomics.com/d3b-bixu/add-strelka2-fields:1.0.0
      - `Mutect2` and all `GATK` tools: pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
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

  ![pipeline flowchart](./docs/kfdrc-somatic-variant-workflow.cwl.png)

requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement
inputs:
  # Required
  reference_fasta: {type: 'File', "sbg:suggestedValue": {class: File, path: 60639014357c3a53540ca7a3,
      name: Homo_sapiens_assembly38.fasta}}
  reference_fai: {type: 'File?', "sbg:suggestedValue": {class: File, path: 60639016357c3a53540ca7af,
      name: Homo_sapiens_assembly38.fasta.fai}}
  reference_dict: {type: 'File?', "sbg:suggestedValue": {class: File, path: 60639019357c3a53540ca7e7,
      name: Homo_sapiens_assembly38.dict}}
  input_tumor_aligned:
    type: File
    secondaryFiles: |
      ${
        var dpath = self.location.replace(self.basename, "")
        if(self.nameext == '.bam'){
          return {"location": dpath+self.nameroot+".bai", "class": "File"}
        }
        else{
          return {"location": dpath+self.basename+".crai", "class": "File"}
        }
      }
    doc: "tumor BAM or CRAM"
  input_tumor_name: string
  input_normal_aligned:
    type: File
    secondaryFiles: |
      ${
        var dpath = self.location.replace(self.basename, "")
        if(self.nameext == '.bam'){
          return {"location": dpath+self.nameroot+".bai", "class": "File"}
        }
        else{
          return {"location": dpath+self.basename+".crai", "class": "File"}
        }
      }
    doc: "normal BAM or CRAM"
  input_normal_name: string
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache",
    "sbg:suggestedValue": {class: File, path: 607713829360f10e3982a425, name: homo_sapiens_vep_93_GRCh38.tar.gz}}
  cfree_chr_len: {type: 'File', doc: "file with chromosome lengths", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051c4, name: hs38_chr.len}}
  cfree_ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC\
      \ to try"}
  cnvkit_annotation_file: {type: 'File', doc: "refFlat.txt file", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051c1, name: refFlat_HG38.txt}}
  hg38_strelka_bed: {type: 'File', doc: "Bgzipped interval bed file. Recommned padding\
      \ 100bp for WXS; Recommend canonical chromosomes for WGS", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051ae, name: hg38_strelka.bed.gz}}
  hg38_strelka_tbi: {type: 'File?', doc: "Tabix index for hg38_strelka_bed", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051aa, name: hg38_strelka.bed.gz.tbi}}
  extra_arg: {type: 'string?', doc: "Add special options to config file, i.e. --max-input-depth\
      \ 10000"}
  strelka2_cores: {type: 'int?', doc: "Adjust number of cores used to run strelka2",
    default: 18}
  mutect2_af_only_gnomad_vcf: {type: 'File', "sbg:suggestedValue": {class: File, path: 5f50018fe4b054958bc8d2e3,
      name: af-only-gnomad.hg38.vcf.gz}}
  mutect2_af_only_gnomad_tbi: {type: 'File?', doc: "Tabix index for mutect2_af_only_gnomad_vcf",
    "sbg:suggestedValue": {class: File, path: 5f50018fe4b054958bc8d2e5, name: af-only-gnomad.hg38.vcf.gz.tbi}}
  mutect2_exac_common_vcf: {type: 'File', "sbg:suggestedValue": {class: File, path: 5f500135e4b0370371c051ad,
      name: small_exac_common_3.hg38.vcf.gz}}
  mutect2_exac_common_tbi: {type: 'File?', doc: "Tabix index for mutect2_exac_common_vcf",
    "sbg:suggestedValue": {class: File, path: 5f500135e4b0370371c051af, name: small_exac_common_3.hg38.vcf.gz.tbi}}
  output_basename: {type: 'string', doc: "String value to use as basename for outputs"}
  wgs_or_wxs: {type: {type: enum, name: wgs_or_wxs, symbols: ["WGS", "WXS"]}, doc: "Select\
      \ if this run is WGS or WXS"}

  # GATK CNV Inputs
  count_panel_of_normals: {type: 'File?', doc: "Path to read-count PoN created by\
      \ the panel workflow. Significantly reduces quality of calling if not provided!",
    'sbg:fileTypes': "HDF5"}
  run_funcotatesegments: {type: 'boolean', doc: "If true, run Funcotator on the called\
      \ copy-ratio segments. This will generate both a simple TSV and a gene list."}
  funcotator_data_sources_tgz: {type: 'File?', doc: "Path to tar.gz containing the\
      \ data sources for Funcotator to create annotations.", 'sbg:fileTypes': "TAR,\
      \ TAR.GZ, TGZ", 'sbg:suggestedValue': {class: File, path: 60e5f8636a504e4e0c6408d8,
      name: funcotator_dataSources.v1.6.20190124s.tar.gz}}
  funcotator_minimum_segment_size: {type: 'int?', doc: "The minimum number of bases\
      \ for a variant to be annotated as a segment. Recommended to be changed only\
      \ for use with FuncotateSegments. If you encounter 'Variant context does not\
      \ represent a copy number segment' error, set this value lower than the length\
      \ of the failed segment."}

  # Optional with One Default
  cfree_threads: {type: 'int?', default: 16, doc: "For ControlFreeC. Recommend 16\
      \ max, as I/O gets saturated after that losing any advantage"}
  cfree_mate_orientation_control: {type: ['null', {type: enum, name: mate_orientation_control,
        symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends),\
      \ RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  cfree_mate_orientation_sample: {type: ['null', {type: enum, name: mate_orientation_sample,
        symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends),\
      \ RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  lancet_ram: {type: 'int?', default: 12, doc: "Adjust in rare circumstances in which\
      \ 12 GB is not enough"}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: [
          "gatk", "grep"]}], default: "gatk", doc: "Choose 'gatk' for SelectVariants\
      \ tool, or 'grep' for grep expression"}
  min_theta2_frac: {type: 'float?', default: 0.01, doc: "Minimum fraction of genome\
      \ with copy umber alterations.  Default is 0.05, recommend 0.01"}
  vardict_cpus: {type: 'int?', default: 9, doc: "Number of CPUs for Vardict to use"}
  vardict_min_vaf: {type: 'float?', default: 0.05, doc: "Min variant allele frequency\
      \ for vardict to consider. Recommend 0.05"}
  vardict_ram: {type: 'int?', default: 18, doc: "GB of RAM to allocate to Vardict\
      \ (hard-capped)"}
  vep_ref_build: {type: 'string?', default: "GRCh38", doc: "Genome ref build used,\
      \ should line up with cache"}

  # Optional with Multiple Defaults (handled in choose_defaults)
  exome_flag: {type: 'string?', doc: "Whether to run in exome mode for callers. Y for\
      \ WXS, N for WGS"}
  lancet_window: {type: 'int?', doc: "Window size for lancet.  Recommend 500 for WGS;\
      \ 600 for exome+"}
  lancet_padding: {type: 'int?', doc: "Recommend 0 if interval file padded already,\
      \ half window size if not. Recommended: 0 for WXS; 300 for WGS"}
  vardict_padding: {type: 'int?', doc: "Padding to add to input intervals, recommend\
      \ 0 if intervals already padded such as in WXS, 150 if not such as in WGS"}
  cnvkit_wgs_mode: {type: 'string?', doc: "for WGS mode, input Y. leave blank for\
      \ WXS/hybrid mode"}
  i_flag: {type: 'string?', doc: "Flag to intersect germline calls on padded regions.\
      \ Use N if you want to skip this or have a WGS run"}

  # Optional
  b_allele: {type: 'File?', doc: "germline calls, needed for BAF.  GATK HC VQSR input\
      \ recommended.  Tool will prefilter for germline and pass if expression given"}
  b_allele_index: {type: 'File?', doc: "Tabix index for b_allele"}
  cfree_coeff_var: {type: 'float?', default: 0.05, doc: "Coefficient of variation\
      \ to set window size.  Default 0.05 recommended"}
  cfree_contamination_adjustment: {type: 'boolean?', doc: "TRUE or FALSE to have ControlFreec\
      \ estimate normal contam"}
  cfree_sex: {type: ['null', {type: enum, name: cfree_sex, symbols: ["XX", "XY"]}],
    doc: "If known, XX for female, XY for male", default: "XX"}
  cnvkit_sex: {type: ['null', {type: enum, name: cnvkit_sex, symbols: ["x", "y"]}],
    doc: "Sex, for simplicity x for female y for male", default: "x"}
  combined_include_expression: {type: 'string?', doc: "Theta2 Purity value: Filter\
      \ expression if vcf has non-PASS combined calls, use as-needed, default set\
      \ for VarDict Java for VarDict", default: FILTER="PASS" && (INFO/STATUS="Germline"
      | INFO/STATUS="StrongSomatic")}
  combined_exclude_expression: {type: 'string?', doc: "Theta2 Purity value: Filter\
      \ expression if vcf has non-PASS combined calls, use as-needed"}
  use_manta_small_indels: {type: 'boolean?', default: false, doc: "Should the program\
      \ use the small indels output from Manta in Strelka2 calling?"}
  learnorientation_memory: {type: 'int?', doc: "GB of memory to allocate to GATK LearnReadOrientationModel;\
      \ defaults to 4 (hard-capped)"}
  getpileup_memory: {type: 'int?', doc: "GB of memory to allocate to GATK GetPileupSummaries;\
      \ defaults to 2 (hard-capped)"}
  filtermutectcalls_memory: {type: 'int?', doc: "GB of memory to allocate to GATK\
      \ FilterMutectCalls; defaults to 4 (hard-capped)"}
  manta_memory: {type: 'int?', doc: "GB of memory to allocate to Manta; defaults to\
      \ 10 (soft-capped)"}
  manta_cores: {type: 'int?', doc: "Number of cores to allocate to Manta; defaults\
      \ to 18"}

  # annotation vars
  genomic_hotspots: {type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing\
      \ hg38 genomic positions corresponding to hotspots", "sbg:suggestedValue": [{
        class: File, path: 607713829360f10e3982a423, name: tert.bed}]}
  protein_snv_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited\
      \ file(s) containing protein names and amino acid positions corresponding to\
      \ hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a426,
        name: protein_snv_cancer_hotspots_v2.tsv}]}
  protein_indel_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited\
      \ file(s) containing protein names and amino acid position ranges corresponding\
      \ to hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a424,
        name: protein_indel_cancer_hotspots_v2.tsv}]}
  bcftools_annot_columns: {type: 'string', doc: "csv string of columns from annotation\
      \ to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_annot_vcf: {type: 'File', doc: "bgzipped annotation vcf file", "sbg:suggestedValue": {
      class: File, path: 5f50018fe4b054958bc8d2e3, name: af-only-gnomad.hg38.vcf.gz}}
  bcftools_annot_vcf_index: {type: 'File', doc: "index of bcftools_annot_vcf", "sbg:suggestedValue": {
      class: File, path: 5f50018fe4b054958bc8d2e5, name: af-only-gnomad.hg38.vcf.gz.tbi}}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to create\
      \ a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to\
      \ add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to\
      \ establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration,\
      \ recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP()\
      \ <= 7\"), \"AF > 0.001\"]"}
  disable_hotspot_annotation: {type: 'boolean?', doc: "Disable Hotspot Annotation\
      \ and skip this task.", default: false}
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}

  # WGS only Fields
  wgs_calling_interval_list: {type: 'File?', doc: "GATK intervals list-style, or bed\
      \ file.  Recommend canocical chromosomes with N regions removed", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051b6, name: wgs_canonical_calling_regions.hg38.bed}}
  lancet_calling_interval_bed: {type: 'File?', doc: "For WGS, highly recommended to\
      \ use CDS bed, and supplement with region calls from strelka2 & mutect2.  Can\
      \ still give calling list as bed if true WGS calling desired instead of exome+",
    "sbg:suggestedValue": {class: File, path: 5f500135e4b0370371c051c0, name: GRCh38.gencode.v31.CDS.merged.bed}}

  # WXS only Fields
  padded_capture_regions: {type: 'File?', doc: "Recommend 100bp pad, for somatic variant"}
  unpadded_capture_regions: {type: 'File?', doc: "Capture regions with NO padding\
      \ for cnv calling"}

outputs:
  ctrlfreec_pval: {type: 'File', outputSource: run_controlfreec/ctrlfreec_pval}
  ctrlfreec_config: {type: 'File', outputSource: run_controlfreec/ctrlfreec_config}
  ctrlfreec_pngs: {type: 'File[]', outputSource: run_controlfreec/ctrlfreec_pngs}
  ctrlfreec_bam_ratio: {type: 'File', outputSource: run_controlfreec/ctrlfreec_bam_ratio}
  ctrlfreec_bam_seg: {type: 'File', outputSource: run_controlfreec/ctrlfreec_bam_seg}
  ctrlfreec_baf: {type: 'File', outputSource: run_controlfreec/ctrlfreec_baf}
  ctrlfreec_info: {type: 'File', outputSource: run_controlfreec/ctrlfreec_info}
  cnvkit_cnr: {type: 'File', outputSource: run_cnvkit/cnvkit_cnr}
  cnvkit_cnn_output: {type: ['null', File], outputSource: run_cnvkit/cnvkit_cnn_output}
  cnvkit_calls: {type: 'File', outputSource: run_cnvkit/cnvkit_calls}
  cnvkit_metrics: {type: 'File', outputSource: run_cnvkit/cnvkit_metrics}
  cnvkit_gainloss: {type: 'File', outputSource: run_cnvkit/cnvkit_gainloss}
  cnvkit_seg: {type: 'File', outputSource: run_cnvkit/cnvkit_seg}
  cnvkit_scatter_plot: {type: 'File', outputSource: run_cnvkit/cnvkit_scatter_plot}
  cnvkit_diagram: {type: 'File', outputSource: run_cnvkit/cnvkit_diagram}
  theta2_calls: {type: 'File?', outputSource: run_theta2_purity/theta2_adjusted_cns}
  theta2_seg: {type: 'File?', outputSource: run_theta2_purity/theta2_adjusted_seg}
  theta2_subclonal_results: {type: ['null', 'File[]'], outputSource: expression_flatten_subclonal_results/output}
  theta2_subclonal_cns: {type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclonal_cns}
  theta2_subclone_seg: {type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclone_seg}
  strelka2_public_outputs: {type: 'File[]', outputSource: run_strelka2/strelka2_public_outputs}
  strelka2_protected_outputs: {type: 'File[]', outputSource: run_strelka2/strelka2_protected_outputs}
  strelka2_prepass_vcf: {type: 'File', outputSource: run_strelka2/strelka2_prepass_vcf}
  manta_pass_vcf: {type: 'File', outputSource: run_manta/manta_pass_vcf}
  manta_prepass_vcf: {type: 'File', outputSource: run_manta/manta_prepass_vcf}
  mutect2_public_outputs: {type: 'File[]', outputSource: run_mutect2/mutect2_public_outputs}
  mutect2_protected_outputs: {type: 'File[]', outputSource: run_mutect2/mutect2_protected_outputs}
  mutect2_prepass_vcf: {type: 'File', outputSource: run_mutect2/mutect2_filtered_vcf}
  vardict_public_outputs: {type: 'File[]', outputSource: run_vardict/vardict_public_outputs}
  vardict_protected_outputs: {type: 'File[]', outputSource: run_vardict/vardict_protected_outputs}
  vardict_prepass_vcf: {type: 'File', outputSource: run_vardict/vardict_prepass_vcf}
  lancet_public_outputs: {type: 'File[]', outputSource: run_lancet/lancet_public_outputs}
  lancet_protected_outputs: {type: 'File[]', outputSource: run_lancet/lancet_protected_outputs}
  lancet_prepass_vcf: {type: 'File', outputSource: run_lancet/lancet_prepass_vcf}
  gatk_cnv_tumor_file_archive: {type: File, outputSource: run_gatk_cnv/tumor_file_archive, doc: "Tar\
      \ archive containing all files generated from the tumor sample aligned reads."}
  gatk_cnv_modeled_segments_tumor: {type: File, outputSource: run_gatk_cnv/modeled_segments_tumor,
    doc: "modelFinal.seg files. These are tab-separated values (TSV) files with a\
      \ SAM-style header containing a read group sample name, a sequence dictionary,\
      \ a row specifying the column headers contained in ModeledSegmentCollection.ModeledSegmentTableColumn,\
      \ and the corresponding entry rows."}
  gatk_cnv_modeled_segments_tumor_plot: {type: File, outputSource: run_gatk_cnv/modeled_segments_tumor_plot,
    doc: "Modeled-segments-plot file. This shows the input denoised copy ratios and/or\
      \ alternate-allele fractions as points, as well as box plots for the available\
      \ posteriors in each segment. The colors of the points alternate with the segmentation.\
      \ Copy ratios are only plotted up to the maximum value specified by the argument\
      \ maximum-copy-ratio. Point sizes can be specified by the arguments point-size-copy-ratio\
      \ and point-size-allele-fraction."}
  gatk_cnv_called_copy_ratio_segments_tumor: {type: File, outputSource: run_gatk_cnv/called_copy_ratio_segments_tumor,
    doc: "Called copy-ratio-segments file. This is a tab-separated values (TSV) file\
      \ with a SAM-style header containing a read group sample name, a sequence dictionary,\
      \ a row specifying the column headers contained in CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn,\
      \ and the corresponding entry rows."}
  gatk_cnv_denoised_tumor_plot: {type: File, outputSource: run_gatk_cnv/denoised_tumor_plot,
    doc: "Denoised-plot file that covers the entire range of the copy ratios"}
  gatk_cnv_normal_file_archive: {type: 'File?', outputSource: run_gatk_cnv/normal_file_archive, doc: "Tar\
      \ archive containing all files generated from the normal sample aligned reads."}
  gatk_cnv_modeled_segments_normal: {type: 'File?', outputSource: run_gatk_cnv/modeled_segments_normal,
    doc: "modelFinal.seg files. These are tab-separated values (TSV) files with a\
      \ SAM-style header containing a read group sample name, a sequence dictionary,\
      \ a row specifying the column headers contained in ModeledSegmentCollection.ModeledSegmentTableColumn,\
      \ and the corresponding entry rows."}
  gatk_cnv_modeled_segments_normal_plot: {type: 'File?', outputSource: run_gatk_cnv/modeled_segments_normal_plot,
    doc: "Modeled-segments-plot file. This shows the input denoised copy ratios and/or\
      \ alternate-allele fractions as points, as well as box plots for the available\
      \ posteriors in each segment. The colors of the points alternate with the segmentation.\
      \ Copy ratios are only plotted up to the maximum value specified by the argument\
      \ maximum-copy-ratio. Point sizes can be specified by the arguments point-size-copy-ratio\
      \ and point-size-allele-fraction."}
  gatk_cnv_called_copy_ratio_segments_normal: {type: 'File?', outputSource: run_gatk_cnv/called_copy_ratio_segments_normal,
    doc: "Called copy-ratio-segments file. This is a tab-separated values (TSV) file\
      \ with a SAM-style header containing a read group sample name, a sequence dictionary,\
      \ a row specifying the column headers contained in CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn,\
      \ and the corresponding entry rows."}
  gatk_cnv_denoised_normal_plot: {type: 'File?', outputSource: run_gatk_cnv/denoised_normal_plot,
    doc: "Denoised-plot file that covers the entire range of the copy ratios"}
  gatk_cnv_funcotated_called_file_tumor: {type: 'File?', outputSource: run_gatk_cnv/funcotated_called_file_tumor,
    doc: "TSV where each row is a segment and the annotations are the covered genes\
      \ and which genes+exon is overlapped by the segment breakpoints."}
  gatk_cnv_funcotated_called_gene_list_file_tumor: {type: 'File?', outputSource: run_gatk_cnv/funcotated_called_gene_list_file_tumor,
    doc: "TSV where each row is a gene and the annotations are the covered genes and\
      \ which genes+exon is overlapped by the segment breakpoints."}

steps:
  choose_defaults:
    run: ../tools/mode_defaults.cwl
    in:
      input_mode: wgs_or_wxs
      exome_flag: exome_flag
      cnvkit_wgs_mode: cnvkit_wgs_mode
      i_flag: i_flag
      lancet_padding: lancet_padding
      lancet_window: lancet_window
      vardict_padding: vardict_padding
    out: [out_exome_flag, out_cnvkit_wgs_mode, out_i_flag, out_lancet_padding, out_lancet_window,
      out_vardict_padding]

  prepare_reference:
    run: ../sub_workflows/prepare_reference.cwl
    in:
      input_fasta: reference_fasta
      input_fai: reference_fai
      input_dict: reference_dict
    out: [indexed_fasta, reference_dict]

  index_b_allele:
    run: ../tools/tabix_index.cwl
    in:
      input_file: b_allele
      input_index: b_allele_index
    out: [output]

  index_strelka_bed:
    run: ../tools/tabix_index.cwl
    in:
      input_file: hg38_strelka_bed
      input_index: hg38_strelka_tbi
    out: [output]

  index_mutect_gnomad:
    run: ../tools/tabix_index.cwl
    in:
      input_file: mutect2_af_only_gnomad_vcf
      input_index: mutect2_af_only_gnomad_tbi
    out: [output]

  index_mutect_exac:
    run: ../tools/tabix_index.cwl
    in:
      input_file: mutect2_exac_common_vcf
      input_index: mutect2_exac_common_tbi
    out: [output]

  index_bcftools_annot_vcf:
    run: ../tools/tabix_index.cwl
    in:
      input_file: bcftools_annot_vcf
      input_index: bcftools_annot_vcf_index
    out: [output]

  select_interval_list:
    run: ../tools/mode_selector.cwl
    in:
      input_mode: wgs_or_wxs
      wgs_input: wgs_calling_interval_list
      wxs_input: padded_capture_regions
    out: [output]

  # WGS only
  python_vardict_interval_split:
    run: ../tools/python_vardict_interval_split.cwl
    doc: "Custom interval list generation for vardict input. Briefly, ~60M bp per\
      \ interval list, 20K bp intervals, lists break on chr and N reginos only"
    in:
      wgs_bed_file: select_interval_list/output
    out: [split_intervals_bed]

  bedtools_intersect_germline:
    run: ../tools/bedtools_intersect.cwl
    in:
      input_vcf: index_b_allele/output
      output_basename: output_basename
      input_bed_file: unpadded_capture_regions
      flag: choose_defaults/out_i_flag
    out: [intersected_vcf]

  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: select_interval_list/output
      reference_dict: prepare_reference/reference_dict
      exome_flag: choose_defaults/out_exome_flag
      scatter_ct:
        valueFrom: ${return 50}
      bands:
        valueFrom: ${return 80000000}
    out: [output]

  gatk_filter_germline:
    run: ../tools/gatk_filter_germline_variant.cwl
    in:
      input_vcf: bedtools_intersect_germline/intersected_vcf
      reference_fasta: prepare_reference/indexed_fasta
      output_basename: output_basename
    out: [filtered_vcf, filtered_pass_vcf]

  samtools_cram2bam_plus_calmd_tumor:
    run: ../tools/samtools_cram2bam_plus_calmd.cwl
    in:
      input_reads: input_tumor_aligned
      threads:
        valueFrom: ${return 16;}
      reference: prepare_reference/indexed_fasta
    out: [bam_file]

  samtools_cram2bam_plus_calmd_normal:
    run: ../tools/samtools_cram2bam_plus_calmd.cwl
    in:
      input_reads: input_normal_aligned
      threads:
        valueFrom: ${return 16;}
      reference: prepare_reference/indexed_fasta
    out: [bam_file]

  select_vardict_bed_interval:
    run: ../tools/mode_list_selector.cwl
    in:
      input_mode: wgs_or_wxs
      wgs_input: python_vardict_interval_split/split_intervals_bed
      wxs_input: gatk_intervallisttools/output
    out: [output]

  run_vardict:
    hints:
    - class: 'sbg:AWSInstanceType'
      value: c5.9xlarge
    run: ../sub_workflows/kfdrc_vardict_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      input_normal_name: input_normal_name
      output_basename: output_basename
      reference_dict: prepare_reference/reference_dict
      bed_invtl_split: select_vardict_bed_interval/output
      padding: choose_defaults/out_vardict_padding
      min_vaf: vardict_min_vaf
      select_vars_mode: select_vars_mode
      cpus: vardict_cpus
      ram: vardict_ram
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_annot_vcf: index_bcftools_annot_vcf/output
      bcftools_public_filter: bcftools_public_filter
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
    out: [vardict_prepass_vcf, vardict_protected_outputs, vardict_public_outputs]

  select_mutect_bed_interval:
    run: ../tools/mode_list_selector.cwl
    in:
      input_mode: wgs_or_wxs
      wgs_input: gatk_intervallisttools/output
      wxs_input: gatk_intervallisttools/output
    out: [output]

  run_mutect2:
    hints:
    - class: 'sbg:AWSInstanceType'
      value: c5.9xlarge
    run: ../sub_workflows/kfdrc_mutect2_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      reference_dict: prepare_reference/reference_dict
      bed_invtl_split: select_mutect_bed_interval/output
      af_only_gnomad_vcf: index_mutect_gnomad/output
      exac_common_vcf: index_mutect_exac/output
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      exome_flag: choose_defaults/out_exome_flag
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      output_basename: output_basename
      learnorientation_memory: learnorientation_memory
      getpileup_memory: getpileup_memory
      filtermutectcalls_memory: filtermutectcalls_memory
      select_vars_mode: select_vars_mode
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_annot_vcf: index_bcftools_annot_vcf/output
      bcftools_public_filter: bcftools_public_filter
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
    out: [mutect2_filtered_stats, mutect2_filtered_vcf, mutect2_protected_outputs,
      mutect2_public_outputs]

  run_strelka2:
    run: ../sub_workflows/kfdrc_strelka2_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      reference_dict: prepare_reference/reference_dict
      hg38_strelka_bed: index_strelka_bed/output
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      manta_small_indels: run_manta/manta_small_indels
      use_manta_small_indels: use_manta_small_indels
      exome_flag: choose_defaults/out_exome_flag
      extra_arg: extra_arg
      strelka2_cores: strelka2_cores
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      output_basename: output_basename
      select_vars_mode: select_vars_mode
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_annot_vcf: index_bcftools_annot_vcf/output
      bcftools_public_filter: bcftools_public_filter
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
    out: [strelka2_prepass_vcf, strelka2_protected_outputs, strelka2_public_outputs]

  bedops_gen_lancet_intervals:
    run: ../tools/preprocess_lancet_intervals.cwl
    in:
      strelka2_vcf:
        source: run_strelka2/strelka2_protected_outputs
        valueFrom: "${var str_vcf = self[1]; return str_vcf;}" # vcf is always item 2 from rename file glob in this context
      mutect2_vcf:
        source: run_mutect2/mutect2_protected_outputs
        valueFrom: "${var mut_vcf = self[1]; return mut_vcf;}"
      ref_bed: lancet_calling_interval_bed
      output_basename: output_basename
    out: [run_bed]

  gatk_intervallisttools_exome_plus:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: bedops_gen_lancet_intervals/run_bed
      reference_dict: prepare_reference/reference_dict
      exome_flag:
        valueFrom: ${return "Y";}
      scatter_ct:
        valueFrom: ${return 50}
      bands:
        valueFrom: ${return 80000000}
    out: [output]

  select_lancet_bed_inteval:
    run: ../tools/mode_list_selector.cwl
    in:
      input_mode: wgs_or_wxs
      wgs_input: gatk_intervallisttools_exome_plus/output
      wxs_input: gatk_intervallisttools/output
    out: [output]

  run_lancet:
    hints:
    - class: 'sbg:AWSInstanceType'
      value: c5.9xlarge
    run: ../sub_workflows/kfdrc_lancet_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      input_normal_name: input_normal_name
      output_basename: output_basename
      select_vars_mode: select_vars_mode
      reference_dict: prepare_reference/reference_dict
      bed_invtl_split: select_lancet_bed_inteval/output
      ram: lancet_ram
      window: choose_defaults/out_lancet_window
      padding: choose_defaults/out_lancet_padding
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_annot_vcf: index_bcftools_annot_vcf/output
      bcftools_public_filter: bcftools_public_filter
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
    out: [lancet_prepass_vcf, lancet_protected_outputs, lancet_public_outputs]

  run_controlfreec:
    run: ../sub_workflows/kfdrc_controlfreec_sub_wf.cwl
    in:
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      input_tumor_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      threads: cfree_threads
      output_basename: output_basename
      ploidy: cfree_ploidy
      mate_orientation_sample: cfree_mate_orientation_sample
      mate_orientation_control: cfree_mate_orientation_control
      capture_regions: unpadded_capture_regions
      indexed_reference_fasta: prepare_reference/indexed_fasta
      reference_fai: reference_fai
      b_allele: gatk_filter_germline/filtered_pass_vcf
      chr_len: cfree_chr_len
      coeff_var: cfree_coeff_var
      contamination_adjustment: cfree_contamination_adjustment
      cfree_sex: cfree_sex
    out: [ctrlfreec_pval, ctrlfreec_config, ctrlfreec_pngs, ctrlfreec_bam_ratio, ctrlfreec_bam_seg,
      ctrlfreec_baf, ctrlfreec_info]

  run_cnvkit:
    run: ../sub_workflows/kfdrc_cnvkit_sub_wf.cwl
    in:
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      tumor_sample_name: input_tumor_name
      input_normal_aligned: samtools_cram2bam_plus_calmd_normal/bam_file
      reference: prepare_reference/indexed_fasta
      normal_sample_name: input_normal_name
      capture_regions: unpadded_capture_regions
      wgs_mode: choose_defaults/out_cnvkit_wgs_mode
      b_allele_vcf: gatk_filter_germline/filtered_pass_vcf
      annotation_file: cnvkit_annotation_file
      output_basename: output_basename
      sex: cnvkit_sex
    out: [cnvkit_cnr, cnvkit_cnn_output, cnvkit_calls, cnvkit_metrics, cnvkit_gainloss,
      cnvkit_seg, cnvkit_scatter_plot, cnvkit_diagram]

  run_theta2_purity:
    run: ../sub_workflows/kfdrc_run_theta2_sub_wf.cwl
    in:
      tumor_cns: run_cnvkit/cnvkit_calls
      reference_cnn: run_cnvkit/cnvkit_cnn_output
      tumor_sample_name: input_tumor_name
      normal_sample_name: input_normal_name
      paired_vcf: run_vardict/vardict_prepass_vcf
      combined_include_expression: combined_include_expression
      combined_exclude_expression: combined_exclude_expression
      min_theta2_frac: min_theta2_frac
      output_basename: output_basename
    out: [theta2_adjusted_cns, theta2_adjusted_seg, theta2_subclonal_results, theta2_subclonal_cns,
      theta2_subclone_seg]

  expression_flatten_subclonal_results:
    run: ../tools/expression_flatten_file_list.cwl
    in:
      input_list: run_theta2_purity/theta2_subclonal_results
    out: [output]

  run_manta:
    run: ../sub_workflows/kfdrc_manta_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      reference_dict: prepare_reference/reference_dict
      hg38_strelka_bed: index_strelka_bed/output
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      vep_cache: vep_cache
      output_basename: output_basename
      manta_memory: manta_memory
      manta_cores: manta_cores
      select_vars_mode: select_vars_mode
    out: [manta_prepass_vcf, manta_pass_vcf, manta_small_indels]

  run_gatk_cnv:
    run: ../sub_workflows/kf_cnv_somatic_pair_wf.cwl
    in:
      input_aligned_reads_tumor: input_tumor_aligned
      input_aligned_reads_normal: input_normal_aligned
      reference_fasta: prepare_reference/indexed_fasta
      reference_dict: prepare_reference/reference_dict
      input_interval_list: select_interval_list/output
      bin_length:
        source: wgs_or_wxs
        valueFrom: |
          $(self == "WGS" ? 1000 : 0)
      padding:
        source: wgs_or_wxs
        valueFrom: |
          $(self == "WGS" ? 250 : 0)
      common_sites: index_b_allele/output
      count_panel_of_normals: count_panel_of_normals
      output_basename: output_basename
      run_funcotatesegments: run_funcotatesegments
      funcotator_data_sources_tgz: funcotator_data_sources_tgz
      funcotator_minimum_segment_size: funcotator_minimum_segment_size
    out: [tumor_file_archive, modeled_segments_tumor, modeled_segments_tumor_plot, called_copy_ratio_segments_tumor, denoised_tumor_plot, normal_file_archive, modeled_segments_normal, modeled_segments_normal_plot, called_copy_ratio_segments_normal, denoised_normal_plot, funcotated_called_file_tumor, funcotated_called_gene_list_file_tumor]

$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: 'sbg:maxNumberOfParallelInstances'
  value: 6
"sbg:license": Apache License 2.0
"sbg:publisher": KFDRC
"sbg:categories":
- BAM
- CNV
- CNVKIT
- CONTROLFREEC
- CRAM
- INDEL
- LANCET
- MAF
- MANTA
- MUTECT2
- SNV
- STRELKA2
- SV
- THETA2
- VARDICT
- VCF
- VEP
sbg:links:
- id: 'https://github.com/kids-first/kf-somatic-workflow/releases/tag/v3.0.1'
  label: github-release
