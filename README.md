Detect cross-contamination and source(s) of the contamination.

Constants used across the code are explained in the R/conta_constants.R.

Tests that can be run with devtools test modules are under the tests folder.

Scripts to run the code are found in scripts folder:
1) Run conta analysis (scripts/conta_run.R):
A single SNPs tsv or pileup file is used as input, along with a dbSNP reference
vcf file to call contamination events. The analysis reports contamination calls,
levels, and plots. It will display cnv metrics and bincounts for Y chromosome
if files are provided.

Input files:
- dbSNP file must contain CAF info field and rsid
- TSV files (two pileup formats are supported, see example inputs under test
folders), must contain chr, pos, and counts for each allele

2) Run source detection (scripts/conta_find_source.R):
This mode requires a set of samples that were already run with the run with
conta analysis. It will use the genotypes for each sample calculated by conta
to find samples that have a likelihood higher than the general likelihood
calculated with the population allele frqeuencies.