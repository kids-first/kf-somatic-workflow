# Kids First DRC Somatic Variant Annotation Workflow
This is a subworkflow that is part of the Kids First DRC Somatic Variant Workflow that can be run as standalone.
Annotation of variant calls helps give context to the possible biological consequences of each variant.

![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

It does the following things as described below:

1. Normalize VCF
1. Strip specified `INFO` and `FORMAT` fields (Only if adding a new annotation that clashes with existing)
1. Annotate with VEP
1. Annotated with an additional vcf - optional, recommend using AF-only WGS gnomad VCF)
1. Soft filter on remarkable variant characteristics
   - KF recommends normal read depth <= 7 and gnomAD AF > 0.001
   - This output will be considered `protected`
1. Annotate with hotspots - KF recommends cancer genome hotspots v2, formatting required and explained below
1. Create MAF output using a modified version of MSKCC's vcf2maf
1. Hard filter on vcf based on user-specified criteria - this output would be considered `public`
1. Create MAF output based on `public` vcf

![annot workflow flowchart](../docs/somatic_annotation_wf.png)

## Workflow Description and KF Recommended Inputs
The additional gnomAD annotation, hotspot annotation, and soft + hard filtering are part of process called "Germline Masking."
The purpose of this is to create outputs that are safe for public consumption by marking and creating a version of outputs deemed a "germline risk" based on specified criteria.
For KF, based on guidance from the Genomic Data Commons (GDC), this means filtering variants with a normal read depth of <= 7 reads and a gnomAD AF > 0.001.
The gnomAD AF filter is pretty intuitive - gnomAD is a database resource of variants and their estimated prevalence in the human population.
Therefore, a variant that is higher than the recommended threshold can be seen as a higher risk of being a common and identifiable variant, and a lower risk for being disease-causing.
The normal depth argument may be less intuitive to some, and an example might help explain its importance:
You've got a somatic variant candidate, and the normal-sample reads all support the reference allele.
However,  if there are only ~5 of those reads, then there's a decent chance that it's a germline het variant and you just didn't get lucky enough to see any alt reads in the normal sample.
A prior understanding that heterozygous germline variant are much more common than somatic variants informs this.

### General workflow inputs:
```yaml
inputs:
  indexed_reference_fasta: {type: 'File', secondaryFiles: [.fai, ^.dict]}
  input_vcf: {type: 'File', secondaryFiles: ['.tbi'], doc: "Input vcf to annotate and soft filter"}
  input_tumor_name: string
  input_normal_name: string
  add_common_fields: {type: 'boolean', doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: false}
  # bcftools strip, if needed
  bcftools_strip_columns: {type: 'string?', doc: "csv string of columns to strip if needed to avoid conflict, i.e INFO/AF"}
  # bcftools annotate if more to do
  bcftools_prefilter_csv: { type: 'string?', doc: "csv of bcftools filter params if you want to prefilter before annotation"}
  bcftools_annot_columns: {type: 'string?', doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF, or syntax to annotate as: INFO/gnomad_3_0_AF:=INFO/AF"}
  bcftools_annot_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "additional bgzipped annotation vcf file"}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration for clues"}
  # VEP-specific
  vep_ram: {type: 'int?', default: 32, doc: "In GB, good for somatic calls, 48+ better for germline"}
  vep_cores: {type: 'int?', default: 16, doc: "Number of cores to use. Good for somatic calls, 32 better for germline"}
  vep_buffer_size: {type: 'int?', doc: "Increase or decrease to balance speed and memory usage. 1000 is good for somatic, 100000 with ram increase better for germline"}
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache"}
  dbnsfp: { type: 'File?', secondaryFiles: [.tbi,^.readme.txt], doc: "VEP-formatted plugin file, index, and readme file containing dbNSFP annotations" }
  dbnsfp_fields: { type: 'string?', doc: "csv string with desired fields to annotate. Use ALL to grab all"}
  merged: { type: 'boolean?', doc: "Set to true if merged cache used", default: true }
  run_cache_existing: { type: boolean, doc: "Run the check_existing flag for cache", default: true }
  run_cache_af: { type: boolean, doc: "Run the allele frequency flags for cache" }
  run_stats: { type: boolean, doc: "Create stats file? Disable for speed", default: false }
  cadd_indels: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD indel annotations" }
  cadd_snvs: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD SNV annotations" }
  # Hotspot Annotation
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task." }
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots" }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots" }
  output_basename: string
  tool_name: string
  # MAF-specific
  retain_info: { type: 'string?', doc: "csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,Hotspot`" }
  retain_fmt: { type: 'string?', doc: "csv string with FORMAT fields that you want to keep" }
  retain_ann: { type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF" }
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}
```

### Recommended reference inputs - all file references can be obtained [here](https://cavatica.sbgenomics.com/u/kfdrc-harmonization/kf-references/)
Secondary files needed for each reference file will be a sub-bullet point
 - `indexed_reference_fasta`: `Homo_sapiens_assembly38.fasta`
   - `Homo_sapiens_assembly38.fasta.fai`
   - `Homo_sapiens_assembly38.dict`
 - `bcftools_strip_columns`: "csv string of columns to strip if needed to avoid conflict, i.e INFO/AF
 - `bcftools_annot_vcf`: `gnomad_3.1.1.vwb_subset.vcf.gz`
   - `gnomad_3.1.1.vwb_subset.vcf.gz.tbi`
 - `bcftools_public_filter`: 'FILTER="PASS"|INFO/HotSpotAllele=1'
 - `bcftools_annot_columns`: "INFO/gnomad_3_1_1_AC:=INFO/AC,INFO/gnomad_3_1_1_AN:=INFO/AN,INFO/gnomad_3_1_1_AF:=INFO/AF,INFO/gnomad_3_1_1_nhomalt:=INFO/nhomalt,INFO/gnomad_3_1_1_AC_popmax:=INFO/AC_popmax,INFO/gnomad_3_1_1_AN_popmax:=INFO/AN_popmax,INFO/gnomad_3_1_1_AF_popmax:=INFO/AF_popmax,INFO/gnomad_3_1_1_nhomalt_popmax:=INFO/nhomalt_popmax,INFO/gnomad_3_1_1_AC_controls_and_biobanks:=INFO/AC_controls_and_biobanks,INFO/gnomad_3_1_1_AN_controls_and_biobanks:=INFO/AN_controls_and_biobanks,INFO/gnomad_3_1_1_AF_controls_and_biobanks:=INFO/AF_controls_and_biobanks,INFO/gnomad_3_1_1_AF_non_cancer:=INFO/AF_non_cancer,INFO/gnomad_3_1_1_primate_ai_score:=INFO/primate_ai_score,INFO/gnomad_3_1_1_splice_ai_consequence:=INFO/splice_ai_consequence" # see https://samtools.github.io/bcftools/bcftools.html#annotate for more tips
 - `gatk_filter_name`: ["NORM_DP_LOW", "GNOMAD_AF_HIGH"]
 - `gatk_filter_expression`: ["vc.getGenotype('insert_normal_sample_name').getDP() <= 7", "gnomad_3_1_1_AF > 0.001"] # NOTE!! Replace `insert_normal_sample_name` with the value you'd use for `input_normal_name`!
 - `vep_cache`: `homo_sapiens_merged_vep_105_GRCh38.tar.gz`
 - `genomic_hotspots`: `tert.bed` # This file has two common TERT promoter gene hot spots
 - `protein_snv_hotspots`: `protein_snv_cancer_hotspots_v2.tsv` # A tsv formatted SNV + MNV subset of https://www.cancerhotspots.org/files/hotspots_v2.xls
 - `protein_indel_hotspots`: `protein_indel_cancer_hotspots_v2.tsv` # A tsv formatted INDEL subset of https://www.cancerhotspots.org/files/hotspots_v2.xls

### Source-specific inputs
For each input, the sub-bullet refers to when to use the suggested input
 - `add_common_fields`
   - Strelka2 calls: `true`
   - All others: `false`
 - `retain_info` # This is fairly subjective, some useful columns unique from each caller to carry over from VCF to MAF
   - Strelka2: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,MQ,MQ0,QSI,HotSpotAllele"
   - Mutect2: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,MBQ,TLOD,HotSpotAllele"
   - Lancet: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,MS,FETS,HotSpotAllele"
   - Vardict: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,MSI,MSILEN,SOR,SSF,HotSpotAllele"
 - `retain_ann` # Similar to above, if run for KF harmonization, recommend the following:
   - Strelka2: "HGVSg"
   - Mutect2: "HGVSg"
   - Lancet: "HGVSg"
   - Vardict: "HGVSg"
 - `bcftools_strip_columns` # if reannotating an old file, especially a KF one that had VEP 93 annotation, recommend the following:
   - "FILTER/GNOMAD_AF_HIGH,FILTER/NORM_DP_LOW,INFO/CSQ,INFO/HotSpotAllele"
 - `bcftools_prefilter_csv` # if annotating a file with calls you want screen for, use this. i.e `FILTER="PASS"`

 - `tool_name`:
   - `Strelka2`: `strelka2_somatic`
   - `Mutect2`: `mutect2_somatic`
   - `Lancet`: `lancet_somatic`
   - `VarDict Java`: `vardict_somatic`
 - `vep_cores`
   - WXS: `16`
   - WGS: `32`
 - `vep_ram`
   - WXS: `32`
   - WGS: `48`
 - `vep_buffer`: `1000`

## Workflow outputs
```yaml
outputs:
  annotated_protected: {type: 'File[]', outputSource: rename_protected/renamed_files}
  annotated_public: {type: 'File[]', outputSource: rename_public/renamed_files}
```
 - `annotated_protected`: `PASS` VCF with annotation pipeline soft `FILTER`-added values, VCF index, and MAF format of VCF
 - `annotated_public_vcf`: Same as `annotated_protected`, hard-filtered to include `PASS` only
