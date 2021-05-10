# Annotation Workflow Spec Upgrade
This workflow is a temp workflow for Ops to patch our current outputs and bring them up to spec.
It will use up to 2 instance per task

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

## Fields that need to be edited
These inputs would be used for every caller. Note: the sample names and base names are EXAMPLES.
Need to edit those for the actual value
```python
input_normal_name = "normal_sample_id"
input_tumor_name = "tumor_sample_id"
output_basename = "task_id"
gatk_filter_expression = ["vc.getGenotype('" + input_normal_name + "').getDP() <= 7", "AF > 0.001"] # replace input_normal_name with normal sample BS ID!!!
```
```yaml
strelka2_pass_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
mutect2_pass_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
lancet_pass_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}
vardict_pass_vcf: {type: File, secondaryFiles: [.tbi], doc: "VEP annotated vcf file."}

```
There are a bunch of other inputs that have sbg suggested files and default values. Don't touch those unless absolutely certain you should!

## Outputs:
```yaml
  annotated_protected_strelka2: {type: 'File[]', outputSource: rename_strelka2_protected/renamed_files}
  annotated_public_strelka2: {type: 'File[]', outputSource: rename_strelka2_public/renamed_files}
  annotated_protected_mutect2: {type: 'File[]', outputSource: rename_mutect2_protected/renamed_files}
  annotated_public_mutect2: {type: 'File[]', outputSource: rename_mutect2_public/renamed_files}
  annotated_protected_lancet: {type: 'File[]', outputSource: rename_lancet_protected/renamed_files}
  annotated_public_lancet: {type: 'File[]', outputSource: rename_lancet_public/renamed_files}
  annotated_protected_vardict: {type: 'File[]', outputSource: rename_vardict_protected/renamed_files}
  annotated_public_vardict: {type: 'File[]', outputSource: rename_vardict_public/renamed_files}
```