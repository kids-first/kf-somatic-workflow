# KFDRC GATK Create CNV Panel Of Normals Workflow

![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

The Kids First Data Resource Center GATK Create CNV Panel of Normals Workflow is
a direct port of the GATK best practices workflow WDL. Use this workflow for
creating a GATK CNV Panel of Normals given a list of normal samples aligned
reads. Supports both WGS and WES/WXS.

## Panel of Normal Recommendations

### Best Case Scenario

Our Panel of Normal standards are the same as those defined in the Broad/GATK
Documentation here:
https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-.

In the best case scenario, you will have 40 samples derived from healthy
tissue. Ideally these samples come from young and healthy individuals (thus
eliminating the chance of a sample containing undiagnosed tumor tissue). Those
40 samples will then be processed in the same way as the tumor sample.
Processing includes all technical variables through sequencing (library
preparation, sequencing platform, etc.). For WGS or WXS that spans the
allosomes, the samples used for the panel of normals should have uniform sex.

### Imperfect Scenarios

Often when working on a project, there will not be 40 normal samples available
to the user. If multiple sequencing approaches/centers are used or you have
allosomal regions of interest, the number of required normals can rapidly
expand. If you are unable to acquire the 40 samples, it is still possible to
create a smaller panel of normals. According to GATK, there is no definitive
rule for the number of samples and even a small panel is better than no panel.

Internally, we have observed that as the panel becomes smaller, more calls are
made and those calls are less accurate. For WGS, we have observed that once
panels become smaller than 15 samples, the calls begin to noticeably
deteriorate. For WXS, calling deteriorates somewhere under 25 samples.

#### Consider Dropping Allosomes

It can be difficult to procure the 80 samples (40 male + 40 female) required to
generate a male and female panel of normals. In case where you have around 40
samples total and perhaps less than 20 male and/or female samples, you might
consider dropping allosomes from the analysis. If your analysis has no regions
of interest on the allosomes, it might be best to drop those chromosomes from
the panel of normals. The mixed-sex panel of normals will be applicable to all
samples in your study at the cost of chrX and chrY calls.

### Run at Own Risk

It is possible to create a panel of normals from a single normal sample.
Internally, we have observed that the calls generated from this scenario are
highly erratic. In the case of WXS, the calls made had little overlap with the
calls made from a larger panel. In the case of WGS, we had instances where the
results looked slightly worse than the 15 sample panel of normals but we also
had instances where the calls in no way resembled those from the larger panel.
Again these results were all better than no panel but were far from what we
would call reliable calls.

#### Potential Workflow

1. Your tumor sample shares an sequencing approach with 40 normals samples (if your analysis includes allosomal regions, these samples must have the same sex as your tumor sample) in the same project.
1. If not, reach out to the sequencing center that provided the tumor sample and request a set of normal samples that share the sequencing approach.
1. If the center cannot provide such samples, search for public or similarly-controlled samples that share the sequencing approach. You will probably have most luck with WGS.

At this point you have exhausted your sample sources. From here you can evaluate where you stand:

- `40 or more`: you meet the GATK recommended minimum. Feel free to proceed.
- `30 to 40`: the results do not meet the GATK minimum but the panel is still rather sizable and worth a run
- `20 to 30`: calls from WXS or low-coverage WGS begin to deteriorate; run at your own risk!
- `10 to 20`: calls from high-coverage WGS begin to deteriorate; run at your own risk!
- `1 to 10`: do not trust calls without extensive follow-up analysis!!!


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

