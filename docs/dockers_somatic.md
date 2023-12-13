# Dockers of kfdrc-somatic-variant-workflow.cwl

TOOL|DOCKER
-|-
aa_classifier.cwl|jluebeck/prepareaa:v0.1203.10
add_strelka2_fields.cwl|pgc-images.sbgenomics.com/d3b-bixu/add-strelka2-fields:1.0.0
amplicon_architect.cwl|jluebeck/prepareaa:v0.1203.10
annotsv.cwl|pgc-images.sbgenomics.com/d3b-bixu/annotsv:3.1.1
awk_chrlen_builder.cwl|ubuntu:22.04
awk_min_seg_length.cwl|ubuntu:20.04
bcbio_filter_vardict_somatic.cwl|pgc-images.sbgenomics.com/d3b-bixu/bcbio_vardict_filter
bcftools_filter_vcf.cwl|pgc-images.sbgenomics.com/d3b-bixu/bvcftools:latest
bcftools_reheader_samples_index.cwl|staphb/bcftools:1.17
bcftools_strip_ann.cwl|pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest
bedtools_intersect.cwl|pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest
clt_pass_file.cwl|None
cns_to_aa_bed.cwl|jluebeck/prepareaa:v0.1203.10
cnvkit_access.cwl|etal/cnvkit:0.9.3
cnvkit_batch.cwl|images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3
cnvkit_batch_only.cwl|images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3
cnvkit_export_theta2.cwl|images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3
cnvkit_import_theta2.cwl|images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3
control-freec-11-6-sbg.cwl|images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1
control_freec_mini_pileup.cwl|images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1
echtvar_anno.cwl|pgc-images.sbgenomics.com/d3b-bixu/echtvar:0.1.9
expression_flatten_file_list.cwl|None
expression_pickvalue_workaround.cwl|None
gatk_Mutect2.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
gatk_bedtointervallist.cwl|broadinstitute/gatk:4.4.0.0
gatk_calculatecontamination.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
gatk_callcopyratiosegments.cwl|broadinstitute/gatk:4.2.4.1
gatk_collectalleliccounts.cwl|broadinstitute/gatk:4.2.4.1
gatk_collectreadcounts.cwl|broadinstitute/gatk:4.2.4.1
gatk_denoisereadcounts.cwl|broadinstitute/gatk:4.2.4.1
gatk_filter_germline_variant.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
gatk_filtermutectcalls.cwl|broadinstitute/gatk:4.1.1.0
gatk_funcotatesegments.cwl|broadinstitute/gatk:4.2.4.1
gatk_gatherpileupsummaries.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
gatk_getpileupsummaries.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
gatk_intervallisttobed.cwl|broadinstitute/gatk:4.4.0.0
gatk_intervallisttools.cwl|broadinstitute/gatk:4.4.0.0
gatk_learnorientationbias.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
gatk_mergemutectstats.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
gatk_mergevcfs.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
gatk_modelsegments.cwl|broadinstitute/gatk:4.2.4.1
gatk_plotdenoisedcopyratios.cwl|broadinstitute/gatk:4.2.4.1
gatk_plotmodeledsegments.cwl|broadinstitute/gatk:4.2.4.1
gatk_preprocessintervals.cwl|broadinstitute/gatk:4.2.4.1
gatk_selectvariants.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
gatk_sortvcf.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
gatk_variant_filter.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
generic_rename_outputs.cwl|None
hotspots_annotation.cwl|quay.io/biocontainers/pysam:0.21.0--py310h41dec4a_1
kf_mskcc_vcf2maf.cwl|pgc-images.sbgenomics.com/d3b-bixu/kf_vcf2maf:v1.0.3
lancet.cwl|pgc-images.sbgenomics.com/d3b-bixu/lancet:1.0.7
manta.cwl|pgc-images.sbgenomics.com/d3b-bixu/manta:1.4.0
normalize_vcf.cwl|pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest
prepare_aa.cwl|jluebeck/prepareaa:v0.1203.10
runtime_validator.cwl|None
samtools_calmd.cwl|pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9
strelka2.cwl|pgc-images.sbgenomics.com/d3b-bixu/strelka
tar.cwl|ubuntu:18.04
theta2_purity.cwl|pgc-images.sbgenomics.com/d3b-bixu/theta2:0.7.1
ubuntu_ratio2seg.cwl|pgc-images.sbgenomics.com/d3b-bixu/python:2.7.13
ubuntu_rename_outputs.cwl|pgc-images.sbgenomics.com/d3b-bixu/ubuntu:18.04
untar_gzip.cwl|ubuntu:20.04
vardictjava.cwl|pgc-images.sbgenomics.com/d3b-bixu/vardict:1.7.0
variant_effect_predictor_105.cwl|ensemblorg/ensembl-vep:release_105.0