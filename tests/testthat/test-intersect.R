# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test tsv and vcf intersect")

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

test_that("Test intersection of pileup with dbSNP to supplementary positions", {

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
