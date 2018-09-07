# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test sample swap")

test_that("conta sample swap detection run", {

  options(warn = -1)

  # Test conta genotype files exists
  expect_true(file.exists(swap_sim_file_1))
  expect_true(file.exists(swap_sim_file_2))
  expect_true(file.exists(swap_sim_file_3))

  # Load the genotype files
  dt1 <- load_conta_file(swap_sim_file_1)
  dt2 <- load_conta_file(swap_sim_file_2)
  dt3 <- load_conta_file(swap_sim_file_3)

  # Test the data tables contain expected data
  expect_equal(nrow(dt1), 10)
  expect_equal(nrow(dt2), 10)
  expect_equal(nrow(dt3), 10)

  # Calculate pairwise concordance
  conc1 <- genotype_concordance(dt1, dt2)
  conc2 <- genotype_concordance(dt1, dt3)
  conc3 <- genotype_concordance(dt2, dt3)

  expect_true(conc1 == 0.9)
  expect_true(conc2 == 0.5)
  expect_true(conc3 == 0.5)
})
