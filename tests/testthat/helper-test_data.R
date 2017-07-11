# Global variables used by tests (helper-* files are loaded by testthat)
data_dir <- "test_data"
out_dir <- "test_output"
dir.create(out_dir, showWarnings = FALSE)

tsvFile = sprintf("%s/test.tsv", data_dir)
pileupFile = sprintf("%s/test.pileup", data_dir)
mafFile = sprintf("%s/test.maf.tsv", out_dir)
mafFile2 = sprintf("%s/test.maf2.tsv", out_dir)
mafFile3 = sprintf("%s/test.maf3.tsv", out_dir)
dbSnpFile = sprintf("%s/test.dbsnp.vcf", data_dir)

errorTsv = sprintf("%s/test.error.tsv", data_dir)
suppTsv = sprintf("%s/test.supp.maf.tsv", data_dir)
wgsTsv = sprintf("%s/test.wgs.maf.tsv", data_dir)
targetedTsv = sprintf("%s/test.targeted.maf.tsv", data_dir)
baseline = sprintf("%s/test.posterior.txt", data_dir)