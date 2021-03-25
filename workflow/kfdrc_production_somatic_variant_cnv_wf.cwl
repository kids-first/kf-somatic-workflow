cwlVersion: v1.0
class: Workflow
id: kfdrc_production_somatic_variant_cnv_wf
label: KFDRC Somatic CNV Workflow
doc: |-
  # KFDRC Somatic CNV Workflow
  Workflow for calling somatic variant calling, copy number variation (CNV), and structural variant
  (SV) calls from a tumor/normal pair of WGS or WXS aligned reads.

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  This workflow takes aligned cram input and performs somatic variant calling using Strelka2, Mutect2,
  Lancet, and VarDict Java, CNV estimation using ControlFreeC and CNVkit, and SV calls using Manta.
  Somatic variant call results are annotated using Variant Effect Predictor, with the Memorial
  Sloane Kettering Cancer Center (MSKCC) vcf2maf wrapper.

  ## Somatic Variant Calling

  [Strelka2](https://github.com/Illumina/strelka) v2.9.3 calls single nucleotide variants (SNV) and insertions/deletions (INDEL).
  [Mutect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php) v4.1.10 from the Broad institute calls SNV, multi-nucleotide variants (MNV, basically equal length substitutions with length > 1) and INDEL.
  [Lancet](https://github.com/nygenome/lancet) v1.0.7 from the New York Genome Center (NYGC) calls SNV, MNV, and INDEL.
  [VarDict Java](https://github.com/AstraZeneca-NGS/VarDictJava) v1.7.0 from AstraZeneca calls SNV, MNV, INDEL and more.
  Each caller has a different approach to variant calling, and together one can glean confident results. Strelka2 is run with default settings, similarly Mutect2 following Broad Best Practices, as of this [workflow](https://github.com/broadinstitute/gatk/blob/4.1.1.0/scripts/mutect2_wdl/mutect2.wdl). Lancet is run in what I'd call an "Exome+" mode, based on the NYGC methods described [here](https://www.biorxiv.org/content/biorxiv/early/2019/04/30/623702.full.pdf). In short, regions from GENCODE gtf with feature annotations `exon`, `UTR`, and start/stop `codon` are used as intervals, as well as regions flanking hits from `strelka2` and `mutect2`. Lastly, VarDict Java run params follow the protocol that the [Blue Collar Bioinformatics](https://bcbio-nextgen.readthedocs.io/en/latest/index.html) uses, with the exception of using a min variant allele frequency (VAF) of 0.05 instead of 0.1, which we find to be relevant for rare cancer variant discovery. We also employ their false positive filtering methods.
  Furthermore, each tool's results, in variant call format (vcf), are filtered on the `PASS` flag, with VarDict Java results additionally filtered for the flag `StrongSomatic`. Their results also include germline hits and other categories by default.
  The pre-`PASS` filtered results can still be obtained from the workflow in the event the user wishes to keep some calls that failed `PASS` criteria.

  ## CNV Estimation

  [ControlFreeC](https://github.com/BoevaLab/FREEC) v11.6 is used for CNV estimation.
  The tool portion of the workflow is a port from the [Seven Bridges Genomics](https://www.sevenbridges.com/) team, with a slight tweak in image outputs.
  Also, the workflow wrapper limits what inputs and outputs are used based on our judgement of utility.
  Outputs include raw ratio calls, copy number calls with p values assigned, b allele frequency data, as well as copy number and b allele frequency plots.
  [CNVkit](https://cnvkit.readthedocs.io/en/stable/) v2.9.3 is a CNV second tool we currently use. [THeTa2](https://github.com/raphael-group/THetA) is used to inform and adjust copy number calls with purity estimations. For both tools, we take advantage of b allele frequency integration for copy number genotype estimation and increased CNV accuracy.

  ## SV Calls

  [Manta](https://github.com/Illumina/manta) v1.4.0 is used to call SVs. Output is also in vcf format, with calls filtered on `PASS`.
  Default settings are used at run time.

  ## Variant Annotation

  [Variant Effect Predictor](https://useast.ensembl.org/info/docs/tools/vep/index.html) release 93, wrapped by [vcf2maf](https://github.com/mskcc/vcf2maf) v1.6.17 is used to annotate somatic variant and SV calls.
  Both the annotated vcf and maf file are made available.

  ## Running WGS or WXS
  This workflow is designed to be able to process either WGS or WXS inputs. This functionality comes from usage of the `wgs_or_wxs`
  input enum. Depending on what is provided for this input, the tool will set the appropriate default values and check that the user
  has provided the correct inputs. For example, if the user sets the input to WGS the lancet_padding value will be defaulted to 300;
  alternatively, if the user sets the input to WXS the lancet_padding value will be defaulted to 0. In either case, the user can
  override the defaults simply by providing their own value for lancet_padding in the inputs.

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

  ### Tips To Run:

  1. For input cram files, be sure to have indexed them beforehand as well.

  1. For ControlFreeC, it is highly recommended that you supply a vcf file with germline calls, GATK Haplotype caller recommended.
  Please also make sure the index for this file is available.
  Also, a range of input ploidy possibilities for the inputs are needed. You can simply use `2`, or put in a range, as an array, like 2, 3, 4.
  For mate orientation, you will need to specify, the drop down and tool doc explains your options.

  1. As a cavatica app, default references for hg38 are already pre-populated, as well as some default settings - i.e., number of threads, coefficient of variation input for ControlFreec, and `PASS` filter tool mode.

  1. What is `select_vars_mode` you ask? On occasion, using GATK's `SelectVariants` tool will fail, so a simple `grep` mode on `PASS` can be used instead.
  Related, `bcftools_filter_vcf` is built in as a convenience in case your b allele frequency file has not been filtered on `PASS`.
  You can use the `include_expression` `Filter="PASS"` to achieve this.

  1. Suggested reference inputs are:

      - `reference_fasta`: [Homo_sapiens_assembly38.fasta](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK
      - `reference_dict`: [Homo_sapiens_assembly38.dict](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK
      - `annotation_file`: [refFlat_HG38.txt](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz) gunzip this file from UCSC.  Needed for gene annotation in `CNVkit`
      - `wgs_calling_interval_list`: [wgs_calling_regions.hg38.interval_list](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1) - need a valid google account, this is a link to the resource bundle from Broad GATK.*To create our 'wgs_canonical_calling_regions.hg38.interval_list', edit this file* by leaving only entries related to chr 1-22, X,Y, and M.M may need to be added.
      - `lancet_calling_interval_bed`: `GRCh38.gencode.v31.CDS.merged.bed`.  As decribed at the beginning, for WGS, it's highly recommended to use CDS bed, and supplement with region calls from strelka2 & mutect2. Our reference was obtained from GENCODE, [release 31](https://www.gencodegenes.org/human/release_31.html) using this gtf file [gencode.v31.primary_assembly.annotation.gtf.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.primary_assembly.annotation.gtf.gz) and parsing features for `UTR`, `start codon`, `stop codon`, and `exon`, then using bedtools sort and merge after converting coordinates into bed format.
      - `padded_capture_regions`: Bed file with exome/targeted/capture regiond.  Recommend 100bp pad, for somatic variant calling
      - `unpadded_capture_regions`: Bed file with exome/targeted/capture regiond.  DO NOT pad, for copy number variation
      - `af_only_gnomad_vcf`: [af-only-gnomad.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/-gatk-best-practices/somatic-hg38) - need a valid google account, this is a link to the best practices google bucket from Broad GATK.
      - `exac_common_vcf`: [small_exac_common_3.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38) - need a valid google account, this is a link to the best practices google bucket from Broad GATK.
      - `hg38_strelka_bed`: [hg38_strelka.bed.gz](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#extended-use-cases) - this link here has the bed-formatted text needed to copy to create this file. You will need to bgzip this file.
       - `vep_cache`: `homo_sapiens_vep_93_GRCh38.tar.gz` from ftp://ftp.ensembl.org/pub/release-93/variation/indexed_vep_cache/ - variant effect predictor cache.
       Current production workflow uses this version, and is compatible with the release used in the vcf2maf tool.
       - `threads`: 16
       - `cfree_chr_len`: hs38_chr.len, this a tsv file with chromosomes and their lengths. Should be limited to canonical chromosomes
        The first column must be chromosomes, optionally the second can be an alternate format of chromosomes.
        Last column must be chromosome length.
        Using the `hg38_strelka_bed`, and removing chrM can be a good source for this.
      - `cfree_coeff_var`: 0.05
      - `cfree_contamination_adjustment`: FALSE

  1. Output files (Note, all vcf files that don't have an explicit index output have index files output as as secondary file.  In other words, they will be captured at the end of the workflow):
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
              - `manta_vep_vcf`: SV call filtered on `PASS`, from manta
              - `manta_vep_tbi`: Index file of above bgzipped vcf
              - `manta_prepass_vcf`: SV results with all `FILTER` categories for manta. Use this file if you believe important variants are being left out when using the algorithm's `PASS` filter.
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
  1. Docker images - the workflow tools will automatically pull them, but as a convenience are listed below:
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
  1. For highly complex samples, some tools have shown themselves to require memory allocation adjustments:
     Manta, GATK LearnReadOrientationModel, GATK GetPileupSummaries, and Vardict. Optional inputs exist to expand the
     memory allocation for these jobs: manta_memory, learnorientation_memory, getpileup_memory, and vardict_ram,
     respectively. For the java tools (Vardict and GATK), these values represent limits on the memory that can
     be used for the respective jobs. Tasks will go to these values and not exceed it. They may succeed or fail,
     but they will not exceed the limit established. The memory allocations for these is hardcapped. The memory
     allocation option for Manta, conversely, is a soft cap. The memory requested will be allocated for the job
     on a particular machine but once the task is running on the machine it may exceed that requested value. For example,
     if Manta's memory allocation is set to 10 GB it will have 10 GB allocated to it at task creation, but, if the
     task ends up running on a machine with more memory available, the task may use it. Setting a value here for Manta
     will not prevent Manta from taking more than that value. The memory usage in Manta is limited by the machine hardware.
     As such the option for Manta memory allocation is described as soft cap. For more information on Manta resource
     usage see their [documentation](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#runtime-hardware-requirements).

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  # Required
  reference_fasta: {type: File }
  reference_fai: { type: 'File?' }
  reference_dict: { type: 'File?' }
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
  vep_cache: {type: File, doc: "tar gzipped cache from ensembl/local converted cache"}
  cfree_chr_len: {type: File, doc: "file with chromosome lengths"}
  cfree_ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try"}
  cnvkit_annotation_file: {type: File, doc: "refFlat.txt file"}
  hg38_strelka_bed: {type: File, doc: "Bgzipped interval bed file. Recommned padding 100bp for WXS; Recommend canonical chromosomes for WGS"}
  hg38_strelka_tbi: {type: 'File?', doc: "Tabix index for hg38_strelka_bed"}
  mutect2_af_only_gnomad_vcf: {type: File}
  mutect2_af_only_gnomad_tbi: {type: 'File?', doc: "Tabix index for mutect2_af_only_gnomad_vcf"}
  mutect2_exac_common_vcf: {type: File}
  mutect2_exac_common_tbi: {type: 'File?', doc: "Tabix index for mutect2_exac_common_vcf"}
  output_basename: {type: string, doc: "String value to use as basename for outputs"}
  wgs_or_wxs: {type: {type: enum, name: sex, symbols: ["WGS", "WXS"] }, doc: "Select if this run is WGS or WXS"}

  # Optional with One Default
  cfree_threads: {type: 'int?', default: 16, doc: "For ControlFreeC. Recommend 16 max, as I/O gets saturated after that losing any advantage"}
  cfree_mate_orientation_control: {type: ['null', {type: enum, name: mate_orientation_control, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  cfree_mate_orientation_sample: {type: ['null', {type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  lancet_ram: {type: 'int?', default: 12, doc: "Adjust in rare circumstances in which 12 GB is not enough"}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], default: "gatk", doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression"}
  min_theta2_frac: {type: 'float?', default: 0.01, doc: "Minimum fraction of genome with copy umber alterations.  Default is 0.05, recommend 0.01"}
  vardict_cpus: {type: 'int?', default: 9, doc: "Number of CPUs for Vardict to use"}
  vardict_min_vaf: {type: 'float?', default: 0.05, doc: "Min variant allele frequency for vardict to consider. Recommend 0.05"}
  vardict_ram: {type: 'int?', default: 18, doc: "GB of RAM to allocate to Vardict (hard-capped)"}
  vep_ref_build: {type: 'string?', default: "GRCh38", doc: "Genome ref build used, should line up with cache"}

  # Optional with Multiple Defaults (handled in choose_defaults)
  exome_flag: {type: string?, doc: "Whether to run in exome mode for callers. Y for WXS, N for WGS"}
  lancet_window: {type: 'int?', doc: "Window size for lancet.  Recommend 500 for WGS; 600 for exome+"}
  lancet_padding: {type: 'int?', doc: "Recommend 0 if interval file padded already, half window size if not. Recommended: 0 for WXS; 300 for WGS"}
  vardict_padding: {type: 'int?', doc: "Padding to add to input intervals, recommend 0 if intervals already padded such as in WXS, 150 if not such as in WGS"}
  cnvkit_wgs_mode: {type: 'string?', doc: "for WGS mode, input Y. leave blank for WXS/hybrid mode"}
  i_flag: {type: 'string?', doc: "Flag to intersect germline calls on padded regions. Use N if you want to skip this or have a WGS run"}

  # Optional
  b_allele: {type: 'File?', doc: "germline calls, needed for BAF.  GATK HC VQSR input recommended.  Tool will prefilter for germline and pass if expression given"}
  b_allele_index: {type: 'File?', doc: "Tabix index for b_allele"}
  cfree_coeff_var: {type: 'float?', default: 0.05, doc: "Coefficient of variation to set window size.  Default 0.05 recommended"}
  cfree_contamination_adjustment: {type: 'boolean?', doc: "TRUE or FALSE to have ControlFreec estimate normal contam"}
  cfree_sex: {type: ['null', {type: enum, name: sex, symbols: ["XX", "XY"] }], doc: "If known, XX for female, XY for male"}
  cnvkit_sex: {type: 'string?', doc: "If known, choices are m,y,male,Male,f,x,female,Female"}
  combined_include_expression: {type: 'string?', doc: "Theta2 Purity value: Filter expression if vcf has non-PASS combined calls, use as-needed, i.e. for VarDict: FILTER=\"PASS\" && (INFO/STATUS=\"Germline\" | INFO/STATUS=\"StrongSomatic\")"}
  combined_exclude_expression: {type: 'string?', doc: "Theta2 Purity value: Filter expression if vcf has non-PASS combined calls, use as-needed"}
  use_manta_small_indels: {type: 'boolean?', default: false, doc: "Should the program use the small indels output from Manta in Strelka2 calling?"}
  learnorientation_memory: {type: 'int?', doc: "GB of memory to allocate to GATK LearnReadOrientationModel; defaults to 4 (hard-capped)"}
  getpileup_memory: {type: 'int?', doc: "GB of memory to allocate to GATK GetPileupSummaries; defaults to 2 (hard-capped)"}
  filtermutectcalls_memory: {type: 'int?', doc: "GB of memory to allocate to GATK FilterMutectCalls; defaults to 4 (hard-capped)"}
  manta_memory: {type: 'int?', doc: "GB of memory to allocate to Manta; defaults to 10 (soft-capped)"}
  manta_cores: {type: 'int?', doc: "Number of cores to allocate to Manta; defaults to 18"}

  # WGS only Fields
  wgs_calling_interval_list: {type: File?, doc: "GATK intervals list-style, or bed file.  Recommend canocical chromosomes with N regions removed"}
  lancet_calling_interval_bed: {type: File?, doc: "For WGS, highly recommended to use CDS bed, and supplement with region calls from strelka2 & mutect2.  Can still give calling list as bed if true WGS calling desired instead of exome+"}

  # WXS only Fields
  padded_capture_regions: {type: 'File?', doc: "Recommend 100bp pad, for somatic variant"}
  unpadded_capture_regions: {type: 'File?', doc: "Capture regions with NO padding for cnv calling"}

outputs:
  ctrlfreec_pval: {type: File, outputSource: run_controlfreec/ctrlfreec_pval}
  ctrlfreec_config: {type: File, outputSource: run_controlfreec/ctrlfreec_config}
  ctrlfreec_pngs: {type: 'File[]', outputSource: run_controlfreec/ctrlfreec_pngs}
  ctrlfreec_bam_ratio: {type: File, outputSource: run_controlfreec/ctrlfreec_bam_ratio}
  ctrlfreec_bam_seg: {type: File, outputSource: run_controlfreec/ctrlfreec_bam_seg}
  ctrlfreec_baf: {type: File, outputSource: run_controlfreec/ctrlfreec_baf}
  ctrlfreec_info: {type: File, outputSource: run_controlfreec/ctrlfreec_info}
  cnvkit_cnr: {type: File, outputSource: run_cnvkit/cnvkit_cnr}
  cnvkit_cnn_output: {type: ['null', File], outputSource: run_cnvkit/cnvkit_cnn_output}
  cnvkit_calls: {type: File, outputSource: run_cnvkit/cnvkit_calls}
  cnvkit_metrics: {type: File, outputSource: run_cnvkit/cnvkit_metrics}
  cnvkit_gainloss: {type: File, outputSource: run_cnvkit/cnvkit_gainloss}
  cnvkit_seg: {type: File, outputSource: run_cnvkit/cnvkit_seg}
  theta2_calls: {type: File?, outputSource: run_theta2_purity/theta2_adjusted_cns}
  theta2_seg: {type: File?, outputSource: run_theta2_purity/theta2_adjusted_seg}
  theta2_subclonal_results: {type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclonal_results}
  theta2_subclonal_cns: {type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclonal_cns}
  theta2_subclone_seg: {type: ['null', 'File[]'], outputSource: run_theta2_purity/theta2_subclone_seg}
  strelka2_vep_vcf: {type: File, outputSource: run_strelka2/strelka2_vep_vcf}
  strelka2_vep_tbi: {type: File, outputSource: run_strelka2/strelka2_vep_tbi}
  strelka2_prepass_vcf: {type: File, outputSource: run_strelka2/strelka2_prepass_vcf}
  strelka2_vep_maf: {type: File, outputSource: run_strelka2/strelka2_vep_maf}
  manta_pass_vcf: {type: File, outputSource: run_manta/manta_pass_vcf}
  manta_prepass_vcf: {type: File, outputSource: run_manta/manta_prepass_vcf}
  mutect2_vep_vcf: {type: File, outputSource: run_mutect2/mutect2_vep_vcf}
  mutect2_vep_tbi: {type: File, outputSource: run_mutect2/mutect2_vep_tbi}
  mutect2_prepass_vcf: {type: File, outputSource: run_mutect2/mutect2_filtered_vcf}
  mutect2_vep_maf: {type: File, outputSource: run_mutect2/mutect2_vep_maf}
  vardict_vep_somatic_only_vcf: {type: File, outputSource: run_vardict/vardict_vep_somatic_only_vcf}
  vardict_vep_somatic_only_tbi: {type: File, outputSource: run_vardict/vardict_vep_somatic_only_tbi}
  vardict_vep_somatic_only_maf: {type: File, outputSource: run_vardict/vardict_vep_somatic_only_maf}
  vardict_prepass_vcf: {type: File, outputSource: run_vardict/vardict_prepass_vcf}
  lancet_vep_vcf: {type: File, outputSource: run_lancet/lancet_vep_vcf}
  lancet_vep_tbi: {type: File, outputSource: run_lancet/lancet_vep_tbi}
  lancet_vep_maf: {type: File, outputSource: run_lancet/lancet_vep_maf}
  lancet_prepass_vcf: {type: File, outputSource: run_lancet/lancet_prepass_vcf}

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
    out: [out_exome_flag,out_cnvkit_wgs_mode,out_i_flag,out_lancet_padding,out_lancet_window,out_vardict_padding]

  prepare_reference:
    run: ../sub_workflows/prepare_reference.cwl
    in:
      input_fasta: reference_fasta
      input_fai: reference_fai
      input_dict: reference_dict
    out: [indexed_fasta,reference_dict]

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
    doc: "Custom interval list generation for vardict input. Briefly, ~60M bp per interval list, 20K bp intervals, lists break on chr and N reginos only"
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
    out:
      [intersected_vcf]

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
    out:
      [filtered_vcf, filtered_pass_vcf]

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
    run: ../tools/mode_selector.cwl
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
    out:
      [vardict_vep_somatic_only_vcf, vardict_vep_somatic_only_tbi, vardict_vep_somatic_only_maf, vardict_prepass_vcf]

  select_mutect_bed_interval:
    run: ../tools/mode_selector.cwl
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
    out:
      [mutect2_filtered_stats, mutect2_filtered_vcf, mutect2_vep_vcf, mutect2_vep_tbi, mutect2_vep_maf]

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
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      output_basename: output_basename
      select_vars_mode: select_vars_mode
    out:
      [strelka2_vep_vcf, strelka2_vep_tbi, strelka2_prepass_vcf, strelka2_vep_maf]

  bedops_gen_lancet_intervals:
    run: ../tools/preprocess_lancet_intervals.cwl
    in:
      strelka2_vcf: run_strelka2/strelka2_vep_vcf
      mutect2_vcf: run_mutect2/mutect2_vep_vcf
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
    run: ../tools/mode_selector.cwl
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
    out:
      [lancet_vep_vcf, lancet_vep_tbi, lancet_vep_maf, lancet_prepass_vcf]

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
    out:
      [ctrlfreec_pval, ctrlfreec_config, ctrlfreec_pngs, ctrlfreec_bam_ratio, ctrlfreec_bam_seg, ctrlfreec_baf, ctrlfreec_info]

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
    out:
      [cnvkit_cnr, cnvkit_cnn_output, cnvkit_calls, cnvkit_metrics, cnvkit_gainloss, cnvkit_seg]

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
    out:
      [theta2_adjusted_cns, theta2_adjusted_seg, theta2_subclonal_results, theta2_subclonal_cns, theta2_subclone_seg]

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
    out:
      [manta_prepass_vcf, manta_pass_vcf, manta_small_indels]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 6
