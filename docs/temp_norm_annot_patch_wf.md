# Annotation Workflow Spec Upgrade
This workflow is a temp workflow for Ops to patch our current outputs and bring them up to spec

## What does it do?
 - Normalize the vcf
 - Run annotation subwf:
    - add standard fields: ONLY FOR STRELKA2
    - VEP annotate
    - Annotate using bcftools to add AF from gnomad
    - Add soft filter, recommend AF < 0.001 and norm sample DP <= 7
    - Annotate hotspots (protected vcf final)
    - Run vcf2maf (protected maf final)
    - Hard filter selecting PASS or HotSpotAllele=1 (public vcf final)
    - Run vcf2maf on that output for public version (public maf final)

## Common inputs
These inputs would be used for every caller. Note: the sample names and base names are EXAMPLES.
Need to edit those for the actual value
```python
input_normal_name = "normal_sample_id"
input_tumor_name = "tumor_sample_id"
output_basename = "task_id"
bcftools_annot_columns = "INFO/AF"
disable_hotspot_annotation = False
gatk_filter_expression = ["vc.getGenotype('" + input_normal_name + "').getDP() <= 7", "AF > 0.001"]
gatk_filter_name = ["NORM_DP_LOW", "GNOMAD_AF_HIGH"]
use_kf_fields = True
bcftools_public_filter = "FILTER=\"PASS\"|INFO/HotSpotAllele=1"

```

## Per tool inputs
### Strelka2
```python
add_common_fields = True
retain_info = "MQ,MQ0,QSI,HotSpotAllele"
tool_name = "strelka2_somatic"
```
### Mutect2
```python
add_common_fields = False
retain_info = "MBQ,TLOD,HotSpotAllele"
tool_name = "mutect2_somatic"
```
### Lancet
```python
add_common_fields = False
retain_info = "MS,FETS,HotSpotAllele"
tool_name = "lancet_somatic"
```
### Vardict
```python
add_common_fields = False
retain_info = "MSI,MSILEN,SOR,SSF,HotSpotAllele"
tool_name = "vardict_somatic"
```

## Outputs:
```yaml
  annotated_protected_vcf: {type: 'File[]', outputSource: rename_pro_vcf/renamed_files}
  annotated_protected_tbi: {type: 'File[]', outputSource: rename_pro_tbi/renamed_files}
  annotated_protected_maf: {type: 'File[]', outputSource: rename_pro_maf/renamed_files}
  annotated_public_vcf: {type: 'File[]', outputSource: rename_pub_vcf/renamed_files}
  annotated_public_tbi: {type: 'File[]', outputSource: rename_pub_tbi/renamed_files}
  annotated_public_maf: {type: 'File[]', outputSource: rename_pub_maf/renamed_files}
```