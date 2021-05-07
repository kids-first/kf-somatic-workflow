# Kids First DRC Strelka2 Variant Calling Workflow
This is a workflow that is the functional equivalent of running only the Strelka2 variant calling portion of the [Kids First DRC Somatic Variant Workflow](../workflow/kfdrc-somatic-variant-workflow.cwl).
Strelka2 is a variant caller that calls single nucleotide variants and small insertions and deletions.

[This workflow](../workflow/kfdrc_production_strelka2_wf.cwl) does the following things as described below:
1. Determine run-time parameters based on whether data is WGS or exome/targeted
1. Create any missing index files
1. Run the [Strelka2 subworkflow](../sub_workflows/kfdrc_strelka2_sub_wf.cwl) as described [here](../docs/kfdrc_strelka2_subworkflow.md)
1. Collect the outputs from the subworkflow
