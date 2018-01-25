# Global variables used by tests (helper-* files are loaded by testthat)
data_dir <- "test_data"
out_dir <- "test_output"
out_dir_wgs <- "test_output/wgs"
out_dir_targeted <- "test_output/targeted"
out_dir_sim <- "test_sim"

# Remove existing test output directories
unlink(out_dir, recursive = TRUE)
unlink(out_dir_sim, recursive = TRUE)

# Create them again
dir.create(out_dir, showWarnings = FALSE)
dir.create(out_dir_wgs, showWarnings = FALSE)
dir.create(out_dir_targeted, showWarnings = FALSE)

tsv_file <- sprintf("%s/test.tsv", data_dir)
pileup_file <- sprintf("%s/test.pileup", data_dir)
male_metrics_file <- sprintf("%s/test.male.bio-metrics.txt", data_dir)
female_metrics_file <- sprintf("%s/test.female.bio-metrics.txt", data_dir)
maf_file <- sprintf("%s/test.maf.tsv", out_dir)
maf_file2 <- sprintf("%s/test.maf2.tsv", out_dir)
maf_file3 <- sprintf("%s/test.maf3.tsv", out_dir)
dbSNP_file <- sprintf("%s/test.dbsnp.vcf", data_dir)
sim_file <- sprintf("%s/test.regular.maf.tsv", data_dir)

error_tsv <- sprintf("%s/test.error.tsv", data_dir)
supp_tsv <- sprintf("%s/test.supp.maf.tsv", data_dir)
wgs_tsv <- sprintf("%s/test.wgs.maf.tsv", data_dir)
targeted_tsv <- sprintf("%s/test.targeted.maf.tsv", data_dir)
baseline <- sprintf("%s/test.posterior.txt", data_dir)

# swap files
swap_test_file_1 <- sprintf("%s/test.sample1_wgs.gt.tsv", data_dir)
swap_test_file_2 <- sprintf("%s/test.sample1_wbc.gt.tsv", data_dir)
swap_test_file_3 <- sprintf("%s/test.sample2_wbc.gt.tsv", data_dir)
