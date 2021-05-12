# Kids First DRC Lancet Variant Calling Workflow
This is a workflow that is the functional equivalent of running only the lancet variant calling portion of the [Kids First DRC Somatic Variant Workflow](../workflow/kfdrc-somatic-variant-workflow.cwl).
Lancet is a variant caller that calls single nucleotide variants, multi-nucleotide variants,  and small insertions and deletions.

[This workflow](../workflow/kfdrc_production_lancet_wf.cwl) does the following things as described below:
1. Determine run-time parameters based on whether data is WGS or exome/targeted
1. Create any missing index files
1. Convert input alignment file to bam and recalculate `MD` tags, if applicable. Recalculating `MD` tags allows for the lancet algorithm to leverage this tag for speed
1. Generate interval lists for scatter jobs
1. If supplied, generate a specialized interval list using a supplied bed file, plus regions from hits in Strelka2 and Mutect2 output
1. Run the [lancet subworkflow](../sub_workflows/kfdrc_lancet_sub_wf.cwl) as described [here](../docs/kfdrc_lancet_sub_wf.md)
1. Collect the outputs from the subworkflow

## Extended Use Case
Lancet uses a graph model to refine INDEL calls, which can be a very time-consuming process. Therefore, lancet is run in what I'd call an "Exome+" mode, based on the NYGC methods described [here](https://www.biorxiv.org/content/biorxiv/early/2019/04/30/623702.full.pdf). In short, regions from GENCODE gtf with feature annotations `exon`, `UTR`, and start/stop `codon` are used as intervals, as well as regions flanking hits from `strelka2` and `mutect2`. This preferred mode of processing WGS input
### Recommended WGS params
 - `exome_flag`: "N"
 - `window`: `600`
 - `padding` `300`
 - `lancet_calling_interval_bed`: `GRCh38.gencode.v31.CDS.merged.bed`

### Recommended WXS params
 - `exome_flag`: "Y"
 - `window`: `600`
 - `padding` `0`
 - `lancet_calling_interval_bed`: `GRCh38.gencode.v31.CDS.merged.bed`
 - `padded_capture_regions`: <bed file with capture regions, recommend padded 100bp>