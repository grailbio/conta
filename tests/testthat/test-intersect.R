# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test conta intersect")

test_that("Test intersection of tsv and dbSNP to generate tsv", {

  conta::intersect_snps(tsv_file, maf_file, dbSNP_file, FALSE)

  expect_true(file.exists(maf_file))

  maf <- read.table(maf_file, header = TRUE, stringsAsFactors = FALSE)
  expect_true(nrow(maf) == 3)
  expect_true(maf$maf[1] == round(0.01078, 4))
  expect_true(maf$maf[2] == round(0.0003994, 4))
  expect_true(maf$maf[3] == round(0.02276, 4))
})

test_that("Test intersection of pileup and dbSNP to generate tsv", {

  conta::intersect_snps(pileup_file, maf_file2, dbSNP_file, FALSE)

  expect_true(file.exists(maf_file2))

  maf <- read.table(maf_file2, header = TRUE, stringsAsFactors = FALSE)
  expect_true(nrow(maf) == 3)
  expect_true(maf$maf[1] == round(0.01078, 4))
  expect_true(maf$maf[2] == round(0.0003994, 4))
  expect_true(maf$maf[3] == round(0.02276, 4))
})

test_that("Test intersection with supplementary positions", {

  conta::intersect_snps(pileup_file, maf_file3, dbSNP_file, TRUE)

  expect_true(file.exists(maf_file3))

  maf <- read.table(maf_file3, header = TRUE, stringsAsFactors = FALSE)
  expect_true(nrow(maf) == 6)
  expect_true(maf$maf[1] == 0)
  expect_true(maf$maf[2] == round(0.01078, 4))
  expect_true(maf$maf[3] == round(0.0003994, 4))
  expect_true(maf$maf[4] == 0)
  expect_true(maf$maf[5] == round(0.02276, 4))
  expect_true(maf$maf[6] == 0)
})

test_that("Test read genome", {

  # Note these are 0-based positions, when intersect is running, it does proper
  # conversion of 1-based to 0-based before calling this function.
  seq <- get_genomic_seq(test_genome, "1", -1, 3, FALSE)
  seq2 <- get_genomic_seq(test_genome, "1", 10, 5, FALSE)
  seq3 <- get_genomic_seq(test_genome, "1", 0, 3, FALSE)
  seq4 <- get_genomic_seq(test_genome, "1", 500, 3, FALSE)
  expect_true(seq == "NNN")
  expect_true(seq2 == "TCACC")
  expect_true(seq3 == "GGG")
  expect_true(seq4 == "NNN")
})

test_that("Test read genome multiple positions", {

  # Note these are 0-based positions, when intersect is running, it does proper
  # conversion of 1-based to 0-based before calling this function.
  seqs <- get_genomic_seqs(test_genome, c("1", "1", "1", "1"),
                           c(-1, 10, 0, 500),
                           c(3, 5, 3, 3),
                           FALSE)

  expect_true(seqs[1] == "NNN")
  expect_true(seqs[2] == "TCACC")
  expect_true(seqs[3] == "GGG")
  expect_true(seqs[4] == "NNN")
})

test_that("Test intersection with context positions", {

  conta::intersect_snps(context_tsv, context_out_tsv, dbSNP_file, FALSE,
                        test_genome)
  expect_true(file.exists(context_out_tsv))

  context <- read.table(context_out_tsv, header = TRUE, stringsAsFactors = FALSE)
  expect_true(nrow(context) == 4)
  expect_true(context$context[1] == "GCC")
  expect_true(context$context[2] == "CCT")
  expect_true(context$context[3] == "GAA")
  expect_true(context$context[4] == "ATG")
})
