# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' @importFrom data.table fread
#' @import conta

context("test sample swap")

test_that("conta sample swap single genotype concordance", {
  options(warn = -1)
  # Test conta genotype files exists
  expect_true(file.exists(swap_sim_file_1))
  expect_true(file.exists(swap_sim_file_2))
  expect_true(file.exists(swap_sim_file_3))
  # Load the genotype files
  dt1 <- data.table::fread(swap_sim_file_1)
  dt2 <- data.table::fread(swap_sim_file_2)
  dt3 <- data.table::fread(swap_sim_file_3)
  # Test the data tables contain expected data
  expect_equal(nrow(dt1), 10)
  expect_equal(nrow(dt2), 10)
  expect_equal(nrow(dt3), 10)
  # Calculate concordance
  conc1 <- conta::genotype_concordance(dt1, dt2)
  conc2 <- conta::genotype_concordance(dt1, dt3)
  conc3 <- conta::genotype_concordance(dt2, dt3)
  expect_true(conc1 == 0.9)
  expect_true(conc2 == 0.5)
  expect_true(conc3 == 0.5)
})

test_that("conta sample swap pairwise genotype concordance", {
  options(warn = -1)
  # Test conta genotype files exists
  expect_true(file.exists(swap_sim_file_1))
  expect_true(file.exists(swap_sim_file_2))
  expect_true(file.exists(swap_sim_file_3))
  # Construct file and label lists
  files <- list(c(swap_sim_file_1, swap_sim_file_2, swap_sim_file_3))
  labels <- list(c("sample1", "sample2", "sample3"))
  # Calculate pairwise concordance, three samples
  concordances <- conta::conta_swap(files, labels)
  # Validate the dimensions of the dataframe
  expect_true(nrow(concordances) == 3)
  expect_true(ncol(concordances) == 3)
  # Validate the expected concordance values
  expect_true(concordances %>% filter(
    Sample1 == "sample1" & Sample2 == "sample2") %>% select(Concordance) == 0.9)
  expect_true(concordances %>% filter(
    Sample1 == "sample1" & Sample2 == "sample3") %>% select(Concordance) == 0.5)
  expect_true(concordances %>% filter(
    Sample1 == "sample2" & Sample2 == "sample3") %>% select(Concordance) == 0.5)
})