## Detect cross-contamination and source

Following scripts are used to run conta toolset:

First install conta library (outside conta folder, run):
R CMD INSTALL --preclean --no-multiarch --with-keep.source conta

### 0) dbSNP file:
Full dbSNP file may be downloaded from:
ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz

### 1) Run contamination analysis (scripts/conta_run.R):
A tsv or pileup file (containing allele counts for each SNP) is used as input,
along with a dbSNP reference vcf file to call contamination events. The analysis
reports contamination calls, levels, and plots. It will display cnv metrics and
bincounts for Y chromosome if files are provided. Input files:
    - dbSNP file must contain CAF info field and rsid.
    - TSV files (two pileup formats are supported, see example inputs under test
folders), must contain chr, pos, and counts for each allele

### 2) Run source detection (scripts/conta_find_source.R):
This mode requires a set of samples that were already run with the run with
conta analysis. It will use the genotypes for each sample calculated by conta
to find samples that have a likelihood higher than the general likelihood
calculated with the population allele frqeuencies.

### 3) Genotype concordance (sample swap) analyses:
Samples that are sequenced from the same genetic donor should have the same
genotypes across SNPs. Conta provides a genotype concordance function to assist
in sample swap analyses. The output of the concordance function is a value
between 0 and 1. Where concordance values close to 1 (above 0.7 in cases where
one of the samples may be contaminated) are considered the same genetic donor.

Expand upon following code to perform pairwise genotype concordance analyses:
```R
conta_gt1 <- load_conta_file("s3:/conta_runs/conta_1/conta_1.gt.tsv")
conta_gt2 <- load_conta_file("s3:/conta_runs/conta_2/conta_2.gt.tsv")
concordance <- genotype_concordance(conta_gt1, conta_gt2)
```

### Output files:
* _<SAMPLE>.conta.tsv_  Contamination quantification main output.
* _<SAMPLE>.bin.lr.png_ Likelihood ratios for each chromosomal regions
* _<SAMPLE>.bin.lr.loh.png_ Likelihood ratios per chromosomal regions (non-LOH)
* _<SAMPLE>.depth.png_    Depth plot pre-filtering of SNPs
* _<SAMPLE>.filtered.depth.png_  Depth plot post-filtering of SNPs
* _<SAMPLE>.error.tsv_    Substitution error model
* _<SAMPLE>.gt.tsv_   Genotype calls
* _<SAMPLE>.gt.loh.tsv_   Genotype calls with LOH regions removed
* _<SAMPLE>.likelihood.png_   Conta maximum likelihood curve
* _<SAMPLE>.log.txt_    log file
* _<SAMPLE>.loh_regions.tsv_    LOH stats for each chromosomal region
* _<SAMPLE>.per_bin.tsv_        Stats for each chromosomal regions
* _<SAMPLE>.per_bin.loh.tsv_    Stats for each chromosomal regions (non-LOH)
* _<SAMPLE>.per_chr.tsv_        Stats for each chromosome
* _<SAMPLE>.per_chr.loh.tsv_    Stats for each chromosome (non-LOH)
* _<SAMPLE>.vfn.cp.png_   Variant frequency (negated) vs. contamination
* _<SAMPLE>.vr.png_     Variant ratios (across sorted locations to visualize LOH)

### Output format:
* **conta_version**   Version of conta that was used
* **conta_call**    Contamination call (tests if avg_log_lr passes a threshold)
* **cf**    Contamination fraction Ignore if conta_call = FALSE
* **sum_log_lr**    Sum of log likelihood ratios across SNPs
* **avg_log_lr**    Average of log likelihood ratios across SNPs
* **snps**    Number of SNPs considered
* **depth**   Mean number of (paired) reads per SNP
* **pos_lr_all**    Fraction of SNPs with positive likelihood ratio
* **pos_lr_x**    Fraction of SNPs with positive likelihood ratio on X chromosome
* **pos_lr_chr_cv**   Coefficient of variation of fraction of SNPs with positive lr
* **y_count**   Fraction of positions on Y chromosome with at least 1 read
* **pregnancy**   Pregnancy call (currently only for male pregnancy if Y chr avail)
* **excluded_regions**    Number of chromosomal regions excluded due to LOH
* **error_rate**    Average substitution error rate per base
* **T>A**    Average specific substitution error rate per base t to a
* **G>A**    Average specific substitution error rate per base g to a
* **C>A**    Average specific substitution error rate per base c to a
* **A>T**    Average specific substitution error rate per base a to t
* **G>T**    Average specific substitution error rate per base g to t
* **C>T**    Average specific substitution error rate per base c to t
* **A>G**    Average specific substitution error rate per base a to g
* **T>G**    Average specific substitution error rate per base t to g
* **C>G**    Average specific substitution error rate per base c to g
* **A>C**    Average specific substitution error rate per base a to c
* **T>C**    Average specific substitution error rate per base t to c
* **G>C**    Average specific substitution error rate per base g to c

### General Guidelines

* Blackswan term is a threshold on the minimum probability a
given event (SNP) may contribute to overall likelihood. Extremely rare events
may get very low probabilities, and this measure prevents one or few
artifactual signals to cause contamination calls. In other terms, blackswan
controls the depth of signal for each SNP.

* Baseline error model (error rate for each loci) may be provided optionally,
otherwise default is to calculate a generic per sample substitution error model.

* To detect contamination with bisulfite converted data, one may use A>T and T>A
SNPs as input (pre-filter dbSNP file), which are unaffected by bisulfite
conversion on CpG contexts. Also allowed are strand specific counts where each
SNP would be counted on a specific strand. See tests for an example.

* Current pregnancy metric can only detect male pregnancy (for female host) by
considering the presence of partial Y chromosome. Y chromosome counts are
provided by biometrics tool. In its absence, this metric will be NA.

