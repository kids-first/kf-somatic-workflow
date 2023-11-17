# Interval Preparation Walkthrough

## Guiding Principle: One Set of Intervals for Every Calling Type

For the sake of post processing analysis it is important that the software for
each analysis in our pipeline are getting as close to the same intervals as
possible. For example, calculating tumor burden is dependent on the total
calling region. If one caller is whole genome and another caller is only coding
regions, they cannot be directly compared. With this caveat in mind, we seek to
provide the same intervals for each calling type: all SNV callers should
receive the same intervals as other SNV callers; all CNV callers should receive
the same intervals as other CNV callers; all SV callers should receive the same
intervals as other SV callers.

This approach might be overly idealistic as each caller has final say on its
calling regions, but, at least from our end, we seek to minimize that source of
difference.

## Getting Started

### Define the Calling Regions

The `calling_regions` are the most foundational inputs for variant calling. In
order to make any calls, the user must define the regions they wish to call.
The regions can be as generic or specific as the user likes. If you're reading
this you likely fall into the group seeking the former. So what are the most
generic calling regions? That depends on how the reads in the alignment file
were created. In general, you're going to want to call all of the regions you
sequenced and aligned. Are your reads from a targeted sequencing experiment or
do they cover the whole genome?

#### WXS Calling Regions

For any targeted sequencing experiment (aka hybridization capture/hybrid
capture/target enrichment), such as WXS, the calling region can be a difficult
concept to tease out. For example, say you are running an exome experiment. It
would make sense that the calling region would just be the regions of the
genome that are defined as exomes. While conceptually correct, this theory
doesn't reflect the actual sequencing experiment.

Targeted sequencing experiment is actually performed using a library of
oligonucleotide baits of ~150 bases that sufficiently cover your target region
(the exome). Therefore, in a WXS experiment the most complete calling region
are the regions covered by the baits. Any WXS experiment should include a BED
file of the baits used to generate the reads. This BED file should be your
input for `calling_regions`.

For more on target vs bait BED files, see [this note](#target-vs-bait-bed-files) from CNVkit.

#### WGS Calling Regions

For whole genome experiments, the calling region is conceptually much more
straightforward. The calling region is the whole genome. Or is it? Unlike with
exons where every region is hand picked to contain easy to sequence and align
reads, whole genome will sequence reads from everywhere: telomeres,
centromeres, repeat regions, etc. These regions, as it turns out, are extremely
difficult if not impossible to align. If you examine the FASTA file used for
alignment, you will notice large swaths of `N`s at the start, end, and
somewhere in the middle of your chromosomes. These are regions that are ignored
during alignment/mapping. These regions traditionally cover telomeres,
centromeres, and other highly repetitive regions; aligning to these regions
fundamentally cannot work with short read NGS; therefore, they are removed from
consideration. Note all FASTA are masked, however, see [this note](#about-ns-in-the-fasta).

How do you get those regions, you ask? Take your reference FASTA and give it to
GATK/Picard ScatterIntervalsByNs. This tool will create an interval list by
splitting a reference at Ns. Once you have those intervals, we also recommend
removing any additional contigs from consideration. WGS calling intervals that
we use were created using the following commands:
```
# - Get the ACGTmers
# - Merge up to 200Ns
gatk ScatterIntervalsByNs
--REFERENCE Homo_sapiens_assembly38.fasta
--OUTPUT calling_regions.interval_list
--OUTPUT_TYPE ACGT
--MAX_TO_MERGE 200

# Then remove everything but chr1:22,X,Y,M
awk -F'\t' '$1 ~ /^@/ || $1 ~ /^chr([12]?[0-9]|[XYM])$/ { print $line }' calling_regions.interval_list > wgs.calling_regions.interval_list
```

Following the above you now have your `calling_regions`. These will work just
fine for SNV and SV calling but they might present problems for CNV calling.
First, `chrM` is not necessary for CNV calling. Second, if your sample is
female, you are also going to want to drop `chrY`. To remove these from
consideration you can use the blacklist...

### Blacklist Regions

Our `calling_regions` have two corresponding blacklist inputs,
`blacklist_regions` and `cnv_blacklist_regions`. These inputs allow the user to
subtract unnecessary or problematic regions from the calling regions. The use
of these inputs is entirely optional but highly recommended for CNV. Above I
describe the alterations that the user should be making to their
`calling_regions` through the use of the `cnv_blacklist_regions` input
(removing `chrM` for WGS, removing `chrY` for females, removing centromeres,
telomeres, etc.).

Most CNV software recommend the removal of centromeres and other difficult
regions; moreover, internal testing has shown that removing high signal regions
can save time and money with SNV and SV tools that do local realignment, such
as VarDict.

### Other Regions Files

There are two more region inputs that affect calling. These files are used for
specific software for specific reasons. I'll keep it brief so we can move
along:

The `coding_sequence_regions` are only relevant for calling WGS data using
Lancet. Presenting Lancet with intervals spanning the whole genome will lead to
extreme runtime and cost. By starting from an limited calling region and
supplementing with calls from Mutect2 and/or Strelka2, Lancet is able to run on
the scale of hours rather than days.

The `chr_len` file is only relevant for calling with Control-FREEC. Control-FREEC
is unique in that its calling regions can only be limited at the chromosome
level. The pipeline will make this file for the user based on the chromosomes
found in the input intervals.

## Interval Processing in the Somatic Workflow

Now that you understand your intervals, you can run the somatic workflow! If you
trust that we're doing everything perfectly, that's great! No need to keep
reading. If not, please allow me to introduce...


### Prepare Regions: Our Interval Processing Engine

At the core of our interval preparation is the
[prepare_regions.cwl](../tools/prepare_regions.cwl) Common Workflow Language (CWL)
workflow. This workflow can perform variable interval manipulation. As it first
starts with a `calling_regions` file. From there what it does depends on user
input:
1. If `calling_padding` provided, pad the calling regions
1. If `supplement_vcfs` provided, concat those records to 1
1. If `blacklist_regions` provided, subtract those regions from 2
1. If `scatter_count` > 0, scatter 3's regions into <scatter_count> files

Any or all of these steps can be skipped. If the user only wishes to apply a
blacklist and scatter, only step 3 onward will be applied to the calling regions.

The workflow also allows the user to dictate fashion in which the intervals are scattered:
- `break_bands_at_multiples_of`: Users can set the bands of the chromosome at
   which there must be a break. For example if the user sets this value to
   1,000,000, any interval that contains a 1,000,000th base of each chromosome
   (chr1:1000000, chr4:100000000, chr22:9000000, etc.) will be split
- `scatter_count`: Users can set the maximum number of scattered intervals to
   make. For example if the user sets this value to 50, IntervalListTools will
   attempt to create 50 or fewer evenly sized files. What's "even" you ask...
- `subdivision_mode`: Users can dictate how the intervals are broken up. This
   can be done with our without breaking the intervals; "even-ness" can be
   based on number of intervals or bases.

The final return of the workflow are:
- processed regions in BED(.GZ) and INTERVALLIST formats
- if scattering, the scattered regions in BED and INTERVALLIST formats

All steps of the workflow are GATK/Picard:
- BedToIntervalList
- IntervalListTools
- IntervalListToBed

### Somatic Workflow Interval Processing Table

| TOOL          | WHOLE GENOME                                                                                                               | WHOLE EXOME                                                                                    |
|---------------|----------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------|
| CNVkit        | accessible regions<br>- blacklist regions                                                                                  | calling regions<br>- blacklist regions                                                         |
| Control-FREEC | chromosomes of calling regions<br>- blacklist regions                                                                      | calling regions<br>- blacklist regions                                                         |
| GATK CNV      | calling regions<br>- cnv blacklist regions<br>From tool: 250 padding                                                       | calling regions<br>- cnv blacklist regions<br>From tool: 250 padding                       |
| Lancet        | coding sequence regions<br>+ Mutect2 VCF<br>+ Strelka2 VCF<br>- calling blacklist<br>+ scattered<br>From tool: 300 padding | calling regions<br>+ 100 padding<br>- blacklist regions<br>From tool: 0 padding                |
| Manta         | calling regions<br>- blacklist regions                                                                                     | calling regions<br>+ 100 padding<br>- blacklist regions                                        |
| Mutect2       | calling regions<br>- blacklist regions<br>+ 80M bands<br>+ scattered                                                       | calling regions<br>+ 100 padding<br>- blacklist regions<br>+ scattered                         |
| Strelka2      | calling regions<br>- blacklist regions                                                                                     | calling regions<br>+ 100 padding<br>- blacklist regions                                        |
| VarDict       | calling regions<br>- blacklist regions<br>+ 20K bands<br>+ scattered<br>From tool: 150 padding                             | calling regions<br>+ 100 padding<br>- blacklist regions<br>+ scattered<br>From tool: 0 padding |

### Whole Genome Sequencing (WGS) Interval Preparation

WGS has three rounds of interval preparation:
1. The [first round](#whole-genome-calling-intervals) creates the intervals for all, except Lancet, SNV and SV callers
1. The [second round](#exome-calling-intervals) creates the intervals for Lancet
1. The [third round](#copy-number-calling-intervals) handles copy number intervals

#### Whole Genome Calling Intervals

The starting point for these intervals is the `calling_regions` [defined above](#wgs-calling-regions).
From there we do the following:
1. start with the `calling_regions`
2. `blacklist_regions` are subtracted from the `calling_regions`
3. the regions from 2 are scattered with bands of 80,000,000 bases
4. the regions from 2 are scattered with bands of 20,000 bases

- Manta uses the unscattered intervals from 2.
- Strelka2 uses the unscattered intervals from 2.
- Mutect2 uses the scattered intervals from 3.
- VarDict uses the scattered intervals from 4.

For more information on the band size please see [this note](#wgs-scatter-band-size-difference).

#### Exome+ Calling Intervals

The second round of interval preparation kicks off after Strelka2 and Mutect2
have finished. The sole purpose of the second round of interval preparation is
to construct a set of calling intervals for Lancet. Lancet is capable of
calling with the true whole genome intervals above; Lancet, however, is
glacially slow when given these intervals. Even with scattering, running Lancet
on the true whole genome intervals can take days. Given this issue, we run
Lancet on what can be described as **EXOME+** intervals. The general idea
behind **EXOME+** is that we supplement coding sequence regions with regions of
interest already identified by Mutect2 or Strelka2. In other words, we're only
using Lancet to check the coding regions and select areas we have reason to
believe variants might exist. Here's what it looks like to make:
1. start with the `coding_sequence_regions`
2. regions from the `strelka2_protected_outputs` and `mutect2_protected_outputs` VCFs are concated to the `coding_sequence_regions`
3. `blacklist_regions` are subtracted from the regions of 2
4. the regions of 3 are scattered without banding

- Lancet uses the scattered intervals from 4.

#### Copy Number Calling Intervals

We ultimately decided to have a separate input for copy number because the copy
number tends to have needs that are not shared by SNV and SV callers.
Additionally, all of our copy number callers have unique approaches to
restricting calling.

##### GATK CNV WGS

The easiest of these to understand is GATK CNV calling. The pipeline GATK has
produced requires that the user provides intervals both at Panel of Normals
(PON) creation as well as at sample processing. For GATK, we provide it both
the `calling_regions` and `cnv_blacklist_regions` as-is.

##### CNVkit

CNVkit handles its calling regions through its [access
command](https://cnvkit.readthedocs.io/en/stable/pipeline.html#access). We run
this command on the reference FASTA, exclude `cnv_blacklist_regions`, and
tolerate streteches of 200 Ns. The resulting file is identical to the
`wgs_calling_regions.interval_list` - `cnv_blacklist_regions`.  The reasoning
behind this removal is because they believe those regions are "inaccessable to
sequencing".

##### Control-FREEC

Control-FREEC provides, ironically, the least control when it comes to
restricting the calling region. The tool provides no way to restrict calling
with intervals. Calling can only be controlled at the chromosome level using a
`chr_len`. Our pipeline will make this file based on the chromosomes present in
the input intervals.

### Whole Exome Sequencing (WXS) Interval Preparation

WXS also has three rounds of interval preparation:
1. The [first round](#snv-sv-bait-capture-regions) creates the intervals for all SNV and SV callers
2. The [second round](#cnv-bait-caputre-regions) creates the capture intervals for Control-FREEC and CNVkit
3. The [third round](#gatk-cnv-wxs) handles the intervals for GATK CNV

#### SNV SV Bait Capture Regionss

As [discussed above](#wxs-calling-regions) exome sequencing starts from the
bait regions. From there, unlike whole genome, we apply a padding of 100 bases
of the intervals. Here's what it looks like to create the `padded capture
regions` for WXS:
1. start with the `calling_regions`
2. pad the calling regions with 100 bases
3. `blacklist_regions` are subtracted from the regions of 2
4. the regions from 2 are scattered without banding

- Manta uses the unscattered bed from 3.
- Strelka2 uses the unscattered bed from 3.
- Lancet uses the scattered beds from 4.
- Mutect2 uses the scattered beds from 4.
- VarDict uses the scattered beds from 4.

#### CNV Bait Capture Regions

A second set of intervals is created exclusively for Control-FREEC and CNVkit.
For WXS, both of these programs have a specific input for capture regions. For
CNVkit, this flag is `--targets`. For Control-FREEC, this flag is
`captureRegions`. In both cases, these programs are expecting the baited
regions [defined above](#wxs-calling-regions). Here's what it looks like to
create the `unpadded caputre regions`:
1. start with the `calling_regions`
2. `blacklist_regions` are subtracted from the regions of 2

- CNVkit uses the unscattered bed from 2.
- Control-FREEC uses the unscattered bed from 2.

#### GATK CNV WXS

The easiest of these to understand is GATK CNV calling. The pipeline GATK has
produced requires that the user provides intervals both at Panel of Normals
(PON) creation as well as at sample processing. For GATK, we provide it both
the `calling_regions` and `cnv_blacklist_regions` as-is.

# Notes

## Target vs Bait BED Files

From the [CNVkit Documentation](https://cnvkit.readthedocs.io/en/master/quickstart.html?#build-a-reference-from-normal-samples-and-infer-tumor-copy-ratios).

> target vs bait BED files: For hybrid capture, the targeted regions (or
> "primary targets") are the genomic regions your capture kit attempts to
> ensure are well covered, e.g. exons of genes of interest. The baited regions
> (or "capture targets") are the genomic regions your kit actually captures,
> usually including about 50bp flanking either side of each target. Give CNVkit
> the bait/capture BED file, not the primary targets.

## About Ns in the FASTA

One of the more interesting things about the N masking assumption (that is:
centromeres, telomeres, and highly repetitive regions are masked in the
reference FASTA), is that it's entirely true for our standard reference FASTA.
GATK's `Homo_sapiens_assembly38.fasta` reference does not N-mask its
centromeres or high signal/repetitive regions found in [ENCODE's
blacklist](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist-README.pdf).

## WGS Scatter Band Size Difference

Band size mostly comes down to computational considerations. The core reason
behind scattering the intervals in the first place is to enable parallel
processing. Typically, when you are dealing with intervals the interval sizes
are not uniform. Having consistent breaks in the chromosome will provide
uniformity, but that can come at the cost of calling. Breaks can come on an
exon or other important calling regions and calling at the ends of intervals is
typically less reliable than in the middle. Finally, there's the consideration
of the maximum interval size your computation can handle. Software performing
memory and/or CPU intensive operations might require small intervals to prevent
runaway resource usage.

Using the above reasoning, we ultimately went with the following for band sizes:
- Mutect2: we chose bands of 80M bases because we wanted to split the intervals
  as little as possible while getting close to 50 scattered files. Without
  bands, our whole genome ends up producing ~25 "balanced" lists. With these 80M
  bands we produce 46 "balanced" lists (a couple of these are very small though).
  The reason we did not go smaller than that is because we noticed we were losing
  calls near the breaks. When we set the bands to 1M we noticed a loss of calls
  that compelled us to go larger and minimize the breaks.
- VarDict: we chose bands of 20K bases because VarDict is extremely
  computationally expensive. Larger bands result cause RAM usage far outside
  what we were willing to allocate to the scattered tasks. 20K was the largest we
  could provide with the default 16 GB of RAM.

## Note from GATK on CNV Intervals

Considerations in genomic intervals are as follows.
- For targeted exomes, the intervals should represent the bait capture or
  target capture regions.
- For whole genomes, either supply regions where coverage is expected across
  samples, e.g. that exclude alternate haplotypes and decoy regions in GRCh38
  or omit the option for references where coverage is expected for the entirety
  of the reference.
- For either type of data, expect to modify the intervals depending on (i)
  extent of masking in the reference used in read mapping and (ii) expectations
  in coverage on allosomal contigs. For example, for mammalian data, expect to
  remove Y chromosome intervals for female samples
