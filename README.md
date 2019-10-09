# KFDRC Somatic Whole Genome Sequence Analysis Workflow

  ![data service logo](https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcS9BnbvIsTkK3QlSGMDvlgu0tZQJ1q4crMvA-S3fcWfIq6y2d2Y)
  
  This is the Kids First Data Resource Center (DRC) Whole Genome Sequencing (WGS) Somatic Workflow, which includes somatic variant calling, copy number variation (CNV), and structural variant (SV) calls. 
  This workflow takes aligned cram input and performs somatic variant calling using Strelka2 and Mutect2, CNV estimation using ControlFreeC, and SV calls using Manta.
  Somatic variant and SV call results are annoated using Variant Effect Predictor, with the Memorial Sloane Kettering Cancer Center (MSKCC) vcf2maf wraper.
  
  ### Somatic Variant Calling:
  [Strelka2](https://github.com/Illumina/strelka) v2.9.3 calls single nucelotide variants (SNVS) and insertions/deletions (INDELS).
  [Mutect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php) v4.1.10 from the Broad institute also calls SNVS and INDELS.
  Each caller has a different approach to variant calling, and together one can glean confident results.  Strelka2 is run with default settings, similarly Mutect2 following Broad Best Practices, as of this [workflow](https://github.com/broadinstitute/gatk/blob/4.1.1.0/scripts/mutect2_wdl/mutect2.wdl).
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
  
  2) For ControlFreeC, it is highly recommended that you supply a vcf file with germline callse, GATK Haplotype caller recommended.
  Please also make sure the index for this file is available.
  Also, a range of input ploidy possibilities for the inputs are needed.  You can simply use `2`, or put in a range, as ann array, like 2, 3, 4.
  
  3) As a cavatica app, default references for hg38 are already pre-populated, as well as some default settings - i.e., number of threads, coefficient of variation input for ControlFreec, and `PASS` filter tool mode.

  4) What is `select_vars_mode` you ask?  On occcasion, using GATK's `SelectVariants` tool will fail, so a simple `grep` mode on `PASS` can be used instead.
  Related, `bcftools_filter_vcf` is built in as a convenience in case youe b allele frequency file has not been filtered on `PASS`.
  You can use the `include_expression` `Filter="PASS"` to achieve this.
  
  5) Suggested inputs are:
```yaml
  reference_fasta: {type: File, doc: "Genome reference fasta file", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe6c', name: 'Homo_sapiens_assembly38.fasta'}}
  reference_dict: {type: File, doc: "Sequence dictionary created using GATK CreateSequenceDictionary from genome fasta file", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe71', name: 'Homo_sapiens_assembly38.dict'}}
  wgs_calling_interval_list: {type: File, doc: "GATK interval calling list", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe6e', name: 'wgs_canonical_calling_regions.hg38.interval_list'}}
  af_only_gnomad_vcf: {type: File, doc: "Broad GATK gnomad reference file", doc: "GATK interval calling list", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe70', name: 'af-only-gnomad.hg38.vcf.gz'}}
  exac_common_vcf: {type: File, doc: "Broad GATK exac reference file", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe72', name: 'small_exac_common_3.hg38.vcf.gz'}}
  hg38_strelka_bed: {type: File, doc: "bgzipped chromosome bed file", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe6d', name: 'hg38_strelka.bed.gz'}}
  threads: 16
  chr_len: {type: File, doc: "file with chromosome lengths", sbg:suggestedValue: {class: 'File', path: '5d9e3424e4b0950cce15fe73', name: 'hs38_chr.len'}}
  coeff_var: 0.05
  contamination_adjustment: FALSE
```
  ### Links/Resources:
  
  The related Github branch for this app is located [here](https://github.com/kids-first/kf-somatic-workflow/tree/mb-publish).
  From here you can pull the project, and modify tool params, inputs, outputs, etc. if you are a developer or a an ambitious researcher.
