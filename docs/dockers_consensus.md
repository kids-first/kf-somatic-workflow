# Dockers of kfdrc_consensus_calling.cwl

TOOL|DOCKER
-|-
add_strelka2_fields.cwl|pgc-images.sbgenomics.com/d3b-bixu/add-strelka2-fields:1.0.0
bcftools_filter_vcf.cwl|pgc-images.sbgenomics.com/d3b-bixu/bvcftools:latest
bcftools_strip_ann.cwl|pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest
consensus_merge.cwl|pgc-images.sbgenomics.com/d3b-bixu/consensus-merge:1.1.0
echtvar_anno.cwl|pgc-images.sbgenomics.com/d3b-bixu/echtvar:0.1.9
gatk_variant_filter.cwl|pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
generic_rename_outputs.cwl|None
hotspots_annotation.cwl|quay.io/biocontainers/pysam:0.21.0--py310h41dec4a_1
kf_mskcc_vcf2maf.cwl|pgc-images.sbgenomics.com/d3b-bixu/kf_vcf2maf:v1.0.3
normalize_vcf.cwl|pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest
prep_mnp_variants.cwl|pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest
variant_effect_predictor_105.cwl|ensemblorg/ensembl-vep:release_105.0
