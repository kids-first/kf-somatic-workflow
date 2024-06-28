# 🏋️‍♂️ Somatic SNV Call Benchmarking Methods and Results
To demonstrate the efficacy and confidence of our caller choice and and consensus calling method, we benchmarked our workflow against a published benchmark dataset.
The dataset was generated in the publication [Establishing community reference samples, data and call sets for benchmarking cancer mutation detection using whole-genome sequencing](https://pubmed.ncbi.nlm.nih.gov/34504347/).
To summarize their work:
Their starting data consists of 21 Replicates distributed among 5 sequencing centers aligned with 3 different aligners (BWA, Bowtie2, NovoAlign). All the aligned BAMs were called with 6 callers (Mutect2, SomaticSniper, VarDict, MuSE, Strelka2, TNscope) resulting into 63 VCFs in total. All the calls were merged together to create a consensus VCF. This consensus file was fed into a deep learning classifier (SomaticSeq, NeuSomatic) that was trained against truth VCFs. These truth VCFs were produced by spiking mutations randomly in the reads. This ML model classify these calls into 4 confidence levels, in which we report the average number of PASS calls given by the classifier in each VCF, as well as REJECT:
 - High Confidence: 60.5 PASS, 0.15 REJECT and high VAF.
 - Medium Confidence: 27.4 PASS, 0.80 REJECT and low VAF not captured in many replicates.
 - Low Confidence: 8.2 PASS, 1.7 REJECT and very low VAF, cannot distinguish from noise.
 - Unclassified: 1.8 PASS and 3.8 REJECT and likely false positive
A helpful [YouTube video](https://www.youtube.com/watch?v=pDsEo0xdHWA) also exists to explain their methodology.
## Comparison to "Gold Standard" Dataset
BAM files relevant to our workflow (BWA-aligned) were called using our standard [somatic workflow](https://github.com/kids-first/kf-somatic-workflow/blob/v4.4.2/workflow/kfdrc-somatic-variant-workflow.cwl) followed up by our [consensus caller](https://github.com/kids-first/kf-somatic-workflow/blob/v4.3.5/workflow/kfdrc_consensus_calling.cwl). Synthetic BAM files did need to have read groups corrected. Originally, their `SM` identifier matched the normal sample used as the starting point for the spike-in experiment. We updated the `SM` field to match the synthetic BAM filename (e.g. `syn_<normal_sample_name>`). The gold standard VCFs provided by the authors had all samples aggregated. To make 1:1 comparisons, we did the following:
1. Gold standard VCFs were downloaded from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/v1.2/ with [SNV](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/v1.2/high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz) and [INDEL](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/v1.2/high-confidence_sINDEL_in_HC_regions_v1.2.vcf.gz) VCFs merged into a single VCF
1. Merged gold standard VCF was then subset to BWA hits and broken up into call sets matching samples (parse `calledSamples` `INFO` field)
1. Next, a custom [benchmark script](https://github.com/kids-first/kfdrc-benchmark/blob/main/scripts/benchmark_vcf_calls.py), which collected the following metrics:
 - True Positive
 - False Positive
 - False Negative
 - F1 Score
 - Precision
 - Accuracy

## Cell Line Consensus Call Method Benchmarking Results
 Results were collated into the following tables, broken down by confidence level aggregations for the cell line data:
### Consensus vs superset (high + med + low + unclassified)

| Sample ID | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| WGS_FD_T_2 | 116706 | 17105 | 5467 | 0.91 | 0.87 | 0.96 |
| WGS_FD_T_3 | 117517 | 16723 | 5763 | 0.91 | 0.88 | 0.95 |
| WGS_NS_T_1 | 110903 | 17673 | 10759 | 0.89 | 0.86 | 0.91 |
| WGS_NS_T_2 | 106410 | 15837 | 10140 | 0.89 | 0.87 | 0.91 |
| WGS_NV_T_2 | 125819 | 12777 | 9314 | 0.92 | 0.91 | 0.93 |
| WGS_NV_T_3 | 125912 | 12079 | 9884 | 0.92 | 0.91 | 0.93 |
| WGS_IL_T_1 | 117910 | 17479 | 4643 | 0.91 | 0.87 | 0.96 |
| WGS_IL_T_2 | 114246 | 15536 | 6339 | 0.91 | 0.88 | 0.95 |

### Consensus vs superset (high + medium + low)

| Sample ID | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| WGS_FD_T_2 | 115883 | 17928 | 5183 | 0.91 | 0.87 | 0.96 |
| WGS_FD_T_3 | 116486 | 17754 | 5424 | 0.91 | 0.87 | 0.96 |
| WGS_NS_T_1 | 110209 | 18367 | 10373 | 0.88 | 0.86 | 0.91 |
| WGS_NS_T_2 | 105843 | 16404 | 9853 | 0.89 | 0.87 | 0.91 |
| WGS_NV_T_2 | 122494 | 16102 | 7134 | 0.91 | 0.88 | 0.94 |
| WGS_NV_T_3 | 122395 | 15596 | 7385 | 0.91 | 0.89 | 0.94 |
| WGS_IL_T_1 | 117510 | 17879 | 4359 | 0.91 | 0.87 | 0.96 |
| WGS_IL_T_2 | 113932 | 15850 | 6046 | 0.91 | 0.88 | 0.95 |

### Consensus vs superset (high + medium)

| Sample ID | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| WGS_FD_T_2 | 111081 | 22730 | 3953 | 0.89 | 0.83 | 0.97 |
| WGS_FD_T_3 | 111394 | 22846 | 4022 | 0.89 | 0.83 | 0.97 |
| WGS_NS_T_1 | 105784 | 22792 | 8613 | 0.87 | 0.82 | 0.92 |
| WGS_NS_T_2 | 102330 | 19917 | 8665 | 0.88 | 0.84 | 0.92 |
| WGS_NV_T_2 | 115284 | 23312 | 3170 | 0.9 | 0.83 | 0.97 |
| WGS_NV_T_3 | 115098 | 22893 | 3311 | 0.9 | 0.83 | 0.97 |
| WGS_IL_T_1 | 113000 | 22389 | 3157 | 0.9 | 0.83 | 0.97 |
| WGS_IL_T_2 | 110223 | 19559 | 4583 | 0.9 | 0.85 | 0.96 |

### Consensus vs superset (high)

| Sample ID | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| WGS_FD_T_2 | 109135 | 24676 | 3644 | 0.89 | 0.82 | 0.97 |
| WGS_FD_T_3 | 109343 | 24897 | 3640 | 0.88 | 0.81 | 0.97 |
| WGS_NS_T_1 | 103771 | 24805 | 8127 | 0.86 | 0.81 | 0.93 |
| WGS_NS_T_2 | 100951 | 21296 | 8303 | 0.87 | 0.83 | 0.92 |
| WGS_NV_T_2 | 112338 | 26258 | 2538 | 0.89 | 0.81 | 0.98 |
| WGS_NV_T_3 | 112229 | 25762 | 2668 | 0.89 | 0.81 | 0.98 |
| WGS_IL_T_2 | 108427 | 21355 | 4150 | 0.89 | 0.84 | 0.96 |
| WGS_IL_T_1 | 110795 | 24594 | 2880 | 0.89 | 0.82 | 0.97 |

## Consensus Synthetic Tumor Benchmark
As above, but this time comparing to the spike-in [synthetic tumor calls](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/DeepLearning_bams/truth_vcf/):

| Sample ID | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| syn_WGS_NS_N_1|  93366 |  3315 |  20872 |  0.89 | 0.97 | 0.82 |
| syn_WGS_NS_N_2|  90758 |  2321 |  23251 |  0.88 | 0.98 | 0.8 |
| syn_WGS_NV_N_2|  99349 |  3318 |  15831 |  0.91 | 0.97 | 0.86 |
| syn_WGS_NV_N_3|  98854 |  3688 |  16146 |  0.91 | 0.96 | 0.86 |
| syn_WGS_IL_N_1|  95173 |  3121 |  19148 |  0.9 | 0.97 | 0.83 |
| syn_WGS_IL_N_2|  93134 |  2730 |  21164 |  0.89 | 0.97 | 0.81 |
| syn_WGS_FD_N_2|  93452 |  3108 |  21614 |  0.88 | 0.97 | 0.81 |
| syn_WGS_FD_N_3|  92335 |  3207 |  22489 |  0.88 | 0.97 | 0.8 |

## Cell Line Per-caller Results
To add further granularity and gauge individual caller quality/contribution to the consensus call method, we also compared results at a caller level, per sample:

### WGS_FD_T_2

| Caller | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Lancet | 103331 | 13511 | 11703 | 0.89 | 0.88 | 0.9 |
| Mutect2 | 78714 | 11239 | 36320 | 0.77 | 0.88 | 0.68 |
| strelka2 | 112891 | 28513 | 2143 | 0.88 | 0.8 | 0.98 |
| Vardict | 110066 | 55041 | 4968 | 0.79 | 0.67 | 0.96 |

### WGS_FD_T_3

| Caller | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Lancet | 103555 | 13415 | 11861 | 0.89 | 0.89 | 0.9 |
| Mutect2 | 78641 | 11283 | 36775 | 0.77 | 0.87 | 0.68 |
| Strelka2 | 113252 | 28603 | 2164 | 0.88 | 0.8 | 0.98 |
| Vardict | 110361 | 55940 | 5055 | 0.78 | 0.66 | 0.96 |

### WGS_IL_T_1

| Caller | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Lancet | 104709 | 12634 | 11448 | 0.9 | 0.89 | 0.9 |
| Mutect2 | 102499 | 18730 | 13658 | 0.86 | 0.85 | 0.88 |
| Strelka2 | 113267 | 26209 | 2890 | 0.89 | 0.81 | 0.98 |
| Vardict | 111830 | 68305 | 4327 | 0.75 | 0.62 | 0.96 |

### WGS_IL_T_2

| Caller | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Lancet | 100939 | 10922 | 13867 | 0.89 | 0.9 | 0.88 |
| Mutect2 | 86329 | 10851 | 28477 | 0.81 | 0.89 | 0.75 |
| Strelka2 | 111108 | 25365 | 3698 | 0.88 | 0.81 | 0.97 |
| Vardict | 110474 | 64339 | 4332 | 0.76 | 0.63 | 0.96 |

### WGS_NS_T_1

| Caller | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Lancet | 95652 | 12003 | 18745 | 0.86 | 0.89 | 0.84 |
| Mutect2 | 87913 | 17439 | 26484 | 0.8 | 0.83 | 0.77 |
| Strelka2 | 112502 | 30198 | 1895 | 0.88 | 0.79 | 0.98 |
| Vardict | 101280 | 73801 | 13117 | 0.7 | 0.58 | 0.89 |

### WGS_NS_T_2

| Caller | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Lancet | 88497 | 9210 | 22498 | 0.85 | 0.91 | 0.8 |
| Mutect2 | 81819 | 13085 | 29176 | 0.79 | 0.86 | 0.74 |
| Strelka2 | 108254 | 27377 | 2741 | 0.88 | 0.8 | 0.98 |
| Vardict | 98589 | 54929 | 12406 | 0.75 | 0.64 | 0.89 |

### WGS_NV_T_2

| Caller | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Lancet | 109614 | 15827 | 8840 | 0.9 | 0.87 | 0.93 |
| Mutect2 | 105573 | 16665 | 12881 | 0.88 | 0.86 | 0.89 |
| Strelka2 | 116451 | 27419 | 2003 | 0.89 | 0.81 | 0.98 |
| Vardict | 113048 | 69223 | 5406 | 0.75 | 0.62 | 0.95 |

### WGS_NV_T_3

| Caller | True Positive | False Positive | False Negative | F1 Score | Precision | Sensitivity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Lancet | 109248 | 15567 | 9161 | 0.9 | 0.88 | 0.92 |
| Mutect2 | 104232 | 15391 | 14177 | 0.88 | 0.87 | 0.88 |
| Strelka2 | 116320 | 27056 | 2089 | 0.89 | 0.81 | 0.98 |
| Vardict | 112944 | 64809 | 5465 | 0.76 | 0.64 | 0.95 |
