# KFDRC GATK Create CNV Panel Of Normals Workflow

![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

The Kids First Data Resource Center GATK Create CNV Panel of Normals Workflow is
a direct port of the GATK best practices workflow WDL. Use this workflow for
creating a GATK CNV Panel of Normals given a list of normal samples aligned
reads. Supports both WGS and WES/WXS.

## Generalized steps:

- prepares a genomic intervals list with PreprocessIntervals
- collects read coverage counts across the preprocessed intervals
- annotate the preprocessed intervals
- creates a CNV PoN with CreateReadCountPanelOfNormals using read coverage
  counts

## Notes from GATK:

- The intervals argument (`input_intervals` and/or `input_interval_list`) is
  required for both WGS and WES workflows and accepts formats compatible with
the GATK -L argument. These intervals will be padded on both sides by the amount
specified by `padding` (default 250) and split into bins of length specified by
`bin_length` (default 1000; specify 0 to skip binning, e.g., for WES).  For WGS,
the intervals should simply cover the autosomal chromosomes (sex chromosomes may
be included, but care should be taken to 1) avoid creating panels of mixed sex,
and 2) denoise case samples only with panels containing only individuals of the
same sex as the case samples).
- Intervals can be blacklisted from coverage collection and all downstream steps
  by using the blacklist intervals argument (`input_exclude_intervals` and/or
`input_exclude_intervals_list`), which accepts formats compatible with the GATK
-XL argument. This may be useful for excluding centromeric regions, etc. from
analysis.  Alternatively, these regions may be manually filtered from the final
callset.

### Other Resources

dockerfiles: https://github.com/d3b-center/bixtools

