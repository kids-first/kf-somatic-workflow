# Kids First DRC Mutect2 Variant Calling Workflow
This is a workflow that is the functional equivalent of running only the Mutect2 variant calling portion of the [Kids First DRC Somatic Variant Workflow](../workflow/kfdrc-somatic-variant-workflow.cwl).
Mutect2 is a variant caller that calls single nucleotide variants, multi-nucleotide variants,  and small insertions and deletions.

[This workflow](../workflow/kfdrc_production_mutect2_wf.cwl) does the following things as described below:
1. Determine run-time parameters based on whether data is WGS or exome/targeted
1. Create any missing index files
1. Generate interval lists for scatter jobs
1. Run the [Mutect2 subworkflow](../sub_workflows/kfdrc_mutect2_sub_wf.cwl) as described [here](../docs/kfdrc_mutect2_sub_wf.md)
1. Collect the outputs from the subworkflow, except for `mutect2_filtered_stats`
