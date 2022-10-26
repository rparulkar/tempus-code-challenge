# tempus-code-challenge

https://github.com/rparulkar/tempus-code-challenge

## Files in the repo

- requirements.txt - contains any non-standard python libraries needed
- annotate-variants.py - tool that queries VEP API to add annotations for variant records
- output.tsv - output from running this tool, with the command below

Run on python3, specifically python 3.9

## Usage
This tool asynchronously scrapes the VEP API to gather annotations for variants with their corresponding HGVS notation.
It takes a source VCF file and output TSV file as its only mandatory arguments. There also exists an optional parameter
to specify assembly. As written, this tool only supports `hg19` and `hg38`. It takes roughly 2 minutes to run with the
given input VCF.

The output contains the following fields:

- Chromosome
- Start Position
- Reference Allele
- Alternate Allele
- Total Depth
- Alternate Allele Depth
- Alternate Allele Frequency
- Minor Allele Frequency (if minor allele matches variant, otherwise None)
- Gene
- ENSEMBL Gene ID
- ENSEMBL Transcript ID
- Canonical Transcript (0 or 1)
- Most Severe Consequence Annotation (from VEP)
- Variant Type (from VEP)

Additionally, the selection criteria for selecting consequences:

1. Select based on `most_severe_consequence` field in annotations
2. If applicable, only output canonical transcripts

## Sample Command

```
annotate_variants.py data/test_vcf_data.txt output.tsv --hg19
```
Assuming the input VCF lives in a folder called `data` in your working directory.
Omission of `--hg19` uses the `hg38` server in the API calls.