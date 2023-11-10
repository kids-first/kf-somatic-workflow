# Echtvar gnomAD Annotation

Our annotation pipeline uses [Echtvar](https://github.com/brentp/echtvar) to
perform gnomAD annotation. Previously, we used BCFtools and observed runtimes of
roughly 30 minutes per file. With Echtvar the same annotations now take under 5
minutes. In order for Echtvar to work, however, it needs to create a custom
reference. The process for creating these custom references is well defined in
the [Echtvar README](https://github.com/brentp/echtvar#encode).

To create our gnomAD annotation file, we used the following files.

## Input VCF

The input VCF was generated from the **genomes** VCFs available from
[gnomAD](https://gnomad.broadinstitute.org/downloads/). From there we selected
only the `PASS` varaints. Additionally, we selected only the following fields:
- `AC`
- `AN`
- `AF`
- `nhomalt`
- `AC_popmax`
- `AN_popmax`
- `AF_popmax`
- `nhomalt_popmax`
- `AC_controls_and_biobanks`
- `AN_controls_and_biobanks`
- `AF_controls_and_biobanks`
- `AF_non_cancer`
- `primate_ai_score`
- `splice_ai_consequence`

Finally, we merged these reads into a singular `gnomad.vwb_subset.vcf.gz`. This
VCF is large (>25GB) but is orders of magnitude smaller than the whole gnomAD
VCF.

In theory, if one had a large enough machine it would be possible to run
Echtvar on the whole genome gnomAD VCFs as detailed in their README but, since
we already had this existing minimal VCF, we went with that. For future
releases of gnomAD we might go with that approach.

## Input JSON

Once you have your VCF, the next step to create the Echtvar reference is to
compose a JSON. The [Echtvar README](https://github.com/brentp/echtvar#configuration-file-for-encode) is
again a great source of information. Below you can see what we ultimately used
for our `gnomad_3_1_1.vwb_subset.echtvar_0_1_9.zip` file. As the versions
change we'll likely update what we use but the general concepts remain true.

```JSON
[
        {"field": "AC", "alias": "gnomad_3_1_1_AC", "missing_value": -2147483648},
        {"field": "AN", "alias": "gnomad_3_1_1_AN", "missing_value": -2147483648},
        {"field": "AF", "alias": "gnomad_3_1_1_AF", "multiplier": 2000000, "missing_value": 2139095041},
        {"field": "nhomalt", "alias": "gnomad_3_1_1_nhomalt", "missing_value": -2147483648},

        {"field": "AC_popmax", "alias": "gnomad_3_1_1_AC_popmax", "missing_value": -2147483648},
        {"field": "AN_popmax", "alias": "gnomad_3_1_1_AN_popmax", "missing_value": -2147483648},
        {"field": "AF_popmax", "alias": "gnomad_3_1_1_AF_popmax", "multiplier": 2000000, "missing_value": 2139095041},
        {"field": "nhomalt_popmax", "alias": "gnomad_3_1_1_nhomalt_popmax", "missing_value": -2147483648},


        {"field": "AC_controls_and_biobanks", "alias": "gnomad_3_1_1_AC_controls_and_biobanks", "missing_value": -2147483648},
        {"field": "AN_controls_and_biobanks", "alias": "gnomad_3_1_1_AN_controls_and_biobanks", "missing_value": -2147483648},
        {"field": "AF_controls_and_biobanks", "alias": "gnomad_3_1_1_AF_controls_and_biobanks", "multiplier": 2000000, "missing_value": 2139095041},

        {"field": "AF_non_cancer", "alias": "gnomad_3_1_1_AF_non_cancer", "multiplier": 2000000, "missing_value": 2139095041},

        {"field": "primate_ai_score", "alias": "gnomad_3_1_1_primate_ai_score", "multiplier": 2000000, "missing_value": 2139095041},
        {"field": "splice_ai_consequence", "alias": "gnomad_3_1_1_splice_ai_consequence", "missing_string": "."},
]
```

1. All `fields` are given a unique `alias`. Here we just prepend all of the
   fields with the release of gnomAD we are using (3.1.1).
1. All integer fields (`AC`, `AN`, `AC_popmax`, etc.) are given a `missing_value` of
   -2147483648. This integer corresponds to `0x80000000` which is the 32 bit
   hexadecimal `missing` value for integers as defined in the [VCF spec](https://samtools.github.io/hts-specs/VCFv4.2.pdf).
1. All float fields (`AF`, `AF_popmax`, `AF_controls_and_biobanks`, etc.) are given a
   `missing_value` of 2139095041. This integer corresponds to `0x7F800001` which
   is the 32 bit hexadecimal `missing` value for floats as defined in the [VCF spec](https://samtools.github.io/hts-specs/VCFv4.2.pdf).
1. All float fields (`AF`, `AF_popmax`, `AF_controls_and_biobanks`, etc.) are given a
   `multiplier` of 2000000. This value is just a sufficiently large value that is
   used to encode/decode the float values. Note that this does end up subtly
   altering the true values albeit very minimally.
