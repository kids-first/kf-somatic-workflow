# KFDRC GATK CNV Somatic Pair Workflow

![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

The Kids First Data Resource Center GATK CNV Somatic Pair Workflow is a direct
port of the GATK best practices workflow WDL. Workflow for running the GATK CNV
pipeline on a matched pair (normal sample optional). Supports both WGS and WES.

## Generalized steps:

- prepares a genomic intervals list with PreprocessIntervals
- collects read coverage counts across the preprocessed intervals
- collects counts of reference versus alternate alleles with
  CollectAllelicCounts
- annotate the preprocessed intervals
- denoises read coverage data against the PoN with DenoiseReadCounts using
  principal component analysis
- plots the results of standardizing and denoising copy ratios against the PoN
- incorporates copy ratio and allelic counts data to group contiguous copy
  ratio and allelic counts segments with ModelSegments using kernel
segmentation and Markov-chain Monte Carlo
- calls amplification, deletion and neutral events for the segmented copy
  ratios
- plots the results of segmentation and estimated allele-specific copy ratios

## Notes from GATK:

- The intervals argument (`input_intervals` and/or `input_interval_list`) is
  required for both WGS and WES workflows and accepts formats compatible with
the GATK -L argument. These intervals will be padded on both sides by the
amount specified by `padding` (default 250) and split into bins of length
specified by `bin_length` (default 1000; specify 0 to skip binning, e.g., for
WES).  For WGS, the intervals should simply cover the autosomal chromosomes
(sex chromosomes may be included, but care should be taken to 1) avoid creating
panels of mixed sex, and 2) denoise case samples only with panels containing
only individuals of the same sex as the case samples).
- Intervals can be blacklisted from coverage collection and all downstream
  steps by using the blacklist intervals argument (`input_exclude_intervals`
and/or `input_exclude_intervals_list`), which accepts formats compatible with
the GATK -XL argument. This may be useful for excluding centromeric regions,
etc. from analysis.  Alternatively, these regions may be manually filtered from
the final callset.
- The sites file (common_sites) should be a Picard or GATK-style interval list.
  This is a list of sites of known variation at which allelic counts will be
collected for use in modeling minor-allele fractions.
- To run `gatk_funcotatesegments.cwl` you must set `run_funcotatesegments` to
  `True` and provide at least `funcotator_data_sources_tgz` and
`funcotator_ref_version`

### Other Resources

dockerfiles: https://github.com/d3b-center/bixtools

