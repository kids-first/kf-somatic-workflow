# Kids First DRC VarDict Java Variant Calling Workflow
This is a workflow that is the functional equivalent of running only the VarDict Java variant calling portion of the [Kids First DRC Somatic Variant Workflow](../workflow/kfdrc-somatic-variant-workflow.cwl).
VarDict Java is a variant caller that calls single nucleotide variants, multi-nucleotide variants,  and small insertions and deletions.

[This workflow](../workflow/kfdrc_production_vardict_wf.cwl) does the following things as described below:
1. Determine run-time parameters based on whether data is WGS or exome/targeted
1. Create any missing index files
1. Generate interval lists for scatter jobs
   - For WGS, will split into 20 kilobase intervals, for a total of 60MB per file for best performance
   - For WXS break into 50 even files
1. If cram input, convert cram to bam
1. Run the [VarDict Java subworkflow](../sub_workflows/kfdrc_vardict_sub_wf.cwl) as described [here](../docs/kfdrc_vardict_sub_wf.md)
1. Collect the outputs from the subworkflow
### Recommended WGS params
 - `vardict_padding`: 150
 - `exome_flag`: "N"
 - `wgs_calling_interval_list`: `wgs_canonical_calling_regions.hg38.bed`
### Recommended WXS params
 - `vardict_padding`: 0
 - `exome_flag`: "Y"
 - `padded_capture_regions`: <bed file with capture regions, recommend padded 100bp>
### Workflow performance tuning
Occasionally, default memory values may not suffice to process the inputs.
Edit where needed to unblock these steps.
All memory units are in whole GB
```yaml
  vardict_cpus: {type: ['null', int], default: 9} # set per task cpu usage
  vardict_ram: {type: ['null', int], default: 18, doc: "In GB"} # set per task memory usage
```
