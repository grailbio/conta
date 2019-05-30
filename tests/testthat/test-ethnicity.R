# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test functions for ethnicity prediction")

test_that("test reading in and parsing vcf", {
  # Load in test vcf and parse. Test vcf file has 10 SNPs per chromosome with
  # setting the removal of sex chromosomes to 'TRUE'.
  parsed_vcf <- parse_1000g_vcf(vcf = "./test_data/test.vcf", remove_sex_chr = TRUE)

  # Expect test vcf's "INFO" column to be split into 10 separate columns,
  # making the parsed_vcf file have 17 columns.
  expect_equal(dim(parsed_vcf)[2], 17)

  # Expect INFO column to separate into specific column names and order
  expect_equal(colnames(parsed_vcf[, 8:17]), c("AC", "AF", "AN", "NS", "DP",
                                               "EAS_AF", "AMR_AF", "AFR_AF",
                                               "EUR_AF", "SAS_AF"))

  # Expect all X and Y chromosome SNPs to be filtered out because the
  # 'parse_1000g_vcf' function when 'remove_sex_chr' is set to TRUE,
  # removes both sex chromosomes
  expect_equal(0, nrow(parsed_vcf %>% dplyr::filter(chr %in% c("X", "Y"))))

  # Because all Y chromosomes are filtered out, expect 220 rows.
  # Rows in test.vcf (240) - X & Y chromosome rows (20) = 220.
  expect_equal(dim(parsed_vcf)[1], 220)

  # There should be no NAs in any column
  expect_equal(FALSE, unique(apply(parsed_vcf, 2, function(x) any(is.na(x)))))

  # Load in test vcf and parse. Test vcf file has 10 SNPs per chromosome with
  # setting the removal of sex chromosomes to 'FALSE'.
  parsed_vcf <- parse_1000g_vcf(vcf = "./test_data/test.vcf", remove_sex_chr = FALSE)

  # Expect all Y chromosome SNPs to be filtered out due to its INFO column
  # having an unexpected ordering of information.
  expect_equal(0, nrow(parsed_vcf %>% dplyr::filter(chr == "Y")))

  # All of chromosome Y was already filtered out due to incorrect ordering of
  # information, however all 10 SNPs on the X chromosome should remain.
  expect_equal(10, nrow(parsed_vcf %>% dplyr::filter(chr == "X")))
})

test_that("test bounding of allele frequencies", {
  # Make fake data with 3 allele allele frequencies
  sim_data  <- c(0, 0.5, 1)

  # Expect that the allele frequency at 0 and 0.5 will increase by 0.1
  # and the allele frequency of 1 will decrease by 0.1.
  new_freqs <- bounded_freqs(sim_data, 0.1)
  expect_equal(new_freqs, c(0.1, 0.6, 0.9))
})

test_that("test prepping of parsed 1000 genome vcf", {
  # Prepare parsed_vcf from the previous test.
  parsed_vcf <- parse_1000g_vcf(vcf = "./test_data/test.vcf")
  prepped_vcf <- prep_1000g_vcf_af(parsed_vcf = parsed_vcf)

  # Expect that there are no population allele frquencies in prepped vcf
  # that equal 0 or 1
  expect_false(any(dplyr::select(prepped_vcf, EAS_AF:SAS_AF) == c(0, 1)))
})

test_that("test receiving snp intersection", {
  # Use prepped_vcf from previous test and take snp intersect with a
  # test.gt.loh.tsv file
  parsed_vcf <- parse_1000g_vcf(vcf = "./test_data/test.vcf")
  prepped_vcf <- prep_1000g_vcf_af(parsed_vcf = parsed_vcf)
  snp_intersect <- get_snp_intersection(gt_loh_file_loc = "./test_data/test.gt.loh.tsv",
                                        prepped_vcf = prepped_vcf)

  # Expect that the intersect between prepped_vcf from the previous test and
  # the test.gt.loh.tsv file has 3 snps
  expect_equal(3, nrow(snp_intersect))
})

test_that("test calculating genotype likelihoods", {
  # Make fake data with 3 genotypes and 3 allele allele frequencies
  sim_data  <- data.frame(genotype = c("0/0", "0/1", "1/1"),
                          af = c(0.1, 0.1, 0.1))

  # Expect that the allele frequency at 0 and 0.5 will increase by 0.1
  # and the allele frequency of 1 will decrease by 0.1.
  gt_likelihoods <- get_gt_likelihoods(genotype = sim_data$genotype,
                                       af = sim_data$af)
  expect_equal(gt_likelihoods, c(0.81, 0.18, 0.01), tolerance = 1e-4)
})

test_that("test calculating snp probabilties", {
  # Use snp_intersect dataframe from previous tests and calculate SNP probabilities
  parsed_vcf <- parse_1000g_vcf(vcf = "./test_data/test.vcf")
  prepped_vcf <- prep_1000g_vcf_af(parsed_vcf = parsed_vcf)
  snp_intersect <- get_snp_intersection(gt_loh_file_loc = "./test_data/test.gt.loh.tsv",
                                        prepped_vcf = prepped_vcf)
  snp_probs <- get_snp_likelihoods(snp_intersect = snp_intersect)

  # Expect that number of rows equals 3 and number of columns equal 23.
  # Columns in snp_intersect (18) + a new column of probabilities per population (5) = 23.
  # Number of rows should equal the number of SNPs in snp_intersect (3).
  expect_equal(dim(snp_probs), c(3, 23))
})

test_that("test ethnicity prediction", {
  # Use initial_snp_likelihoods dataframe from previous tests and calculate SNP probabilities
  parsed_vcf <- parse_1000g_vcf(vcf = "./test_data/test.vcf")
  prepped_vcf <- prep_1000g_vcf_af(parsed_vcf = parsed_vcf)
  snp_intersect <- get_snp_intersection(gt_loh_file_loc =
                                          "./test_data/test.gt.loh.tsv",
                                        prepped_vcf = prepped_vcf)
  initial_snp_likelihoods <- get_snp_likelihoods(snp_intersect)
  ranked_predictions <- calculate_ranked_predictions(initial_snp_likelihoods)

  # Test genotype file was a subset of a self-reported White, Non-Hispanic
  # patient's genotype file (medrio_id = 785), Top ethnicity predictions per
  # selected test chromosomes should be european.
  expect_equal("European",
               unique(ranked_predictions[which(ranked_predictions$pred_rank == 1), ]$ethnicity))
})
