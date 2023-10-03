# Kids First DRC Somatic Variant Annotation Workflow
This is a subworkflow that is part of the Kids First DRC Somatic Variant Workflow that can be run as standalone.
Annotation of variant calls helps give context to the possible biological consequences of each variant.

![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

It does the following things as described below:

1. Normalize VCF
1. Strip specified `INFO` and `FORMAT` fields (Only if adding a new annotation that clashes with existing)
1. Annotate with VEP - can be skipped if VEP run previously and downstream tools are to be repeated
1. Annotated with an additional vcf - optional, recommend using a gnomAD VCF with at least AF
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

### Recommended reference inputs - all file references can be obtained [here](https://cavatica.sbgenomics.com/u/kfdrc-harmonization/kf-references/)
Secondary files needed for each reference file will be a sub-bullet point
 - `indexed_reference_fasta`: `Homo_sapiens_assembly38.fasta`
   - `Homo_sapiens_assembly38.fasta.fai`
   - `Homo_sapiens_assembly38.dict`
 - `echtvar_anno_zips`: `gnomad_3_1_1.vwb_subset.echtvar_0_1_9.zip`
 - `bcftools_strip_columns`: "csv string of columns to strip if needed to avoid conflict, i.e INFO/AF
 - `bcftools_public_filter`: 'FILTER="PASS"|INFO/HotSpotAllele=1'
 - `gatk_filter_name`: ["NORM_DP_LOW", "GNOMAD_AF_HIGH"]
 - `gatk_filter_expression`: ["vc.getGenotype('insert_normal_sample_name').getDP() <= 7", "gnomad_3_1_1_AF != '.' && gnomad_3_1_1_AF > 0.001"] # NOTE!! Replace `insert_normal_sample_name` with the value you'd use for `input_normal_name`! # NOTE!! For the AF filtration to pass, dot values must first be excluded! If they are not, GATK will error trying to convert those values!
 - `vep_cache`: `homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz`
 - `genomic_hotspots`: `tert.bed` # This file has two common TERT promoter gene hot spots
 - `protein_snv_hotspots`: `protein_snv_cancer_hotspots_v2.ENS105_liftover.tsv` # A tsv formatted SNV + MNV subset of https://www.cancerhotspots.org/files/hotspots_v2.xls
 - `protein_indel_hotspots`: `protein_indel_cancer_hotspots_v2.ENS105_liftover.tsv` # A tsv formatted INDEL subset of https://www.cancerhotspots.org/files/hotspots_v2.xls
 - `custom_enst`: `kf_isoform_override.tsv` # As of VEP 104, several genes have had their canonical transcripts redefined. While the VCF will have all possible isoforms, this affects maf file output and may results in representative protein changes that defy historical expectations

### Source-specific inputs
For each input, the sub-bullet refers to when to use the suggested input
 - `add_common_fields`
   - Strelka2 calls: `true`, *exception if already run previously and other downstream tools are being run*
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
 - `bcftools_strip_columns` # if reannotating an old file:
   - "FILTER/GNOMAD_AF_HIGH,FILTER/NORM_DP_LOW,INFO/CSQ,INFO/HotSpotAllele" # recommended if re-annotating from an older VEP cache
   - "FILTER/GNOMAD_AF_HIGH,FILTER/NORM_DP_LOW,INFO/HotSpotAllele" # recommended if repeating hot spot an d want to keep VEP
 - `bcftools_prefilter_csv` # if annotating a file with calls you want screen for, use this. i.e `FILTER="PASS"`
 - `disable_vep_annotation` # set to `True` if existing VEP annotation of file is ok
 - `tool_name`:
   - `Strelka2`: `strelka2_somatic`
   - `Mutect2`: `mutect2_somatic`
   - `Lancet`: `lancet_somatic`
   - `VarDict Java`: `vardict_somatic`
 - `vep_cores`: `16`
 - `vep_ram`: `32`
 - `vep_buffer`: `5000`

## Workflow outputs
 - `annotated_protected`: `PASS` VCF with annotation pipeline soft `FILTER`-added values, VCF index, and MAF format of VCF
 - `annotated_public_vcf`: Same as `annotated_protected`, hard-filtered to include `PASS` only
