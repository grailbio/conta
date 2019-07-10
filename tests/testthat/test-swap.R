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
  expect_true(conc1$Concordance == 0.9)
  expect_true(conc2$Concordance == 0.5)
  expect_true(conc3$Concordance == 0.5)
  expect_true(conc1$No_Change == 9)
  expect_true(conc2$No_Change == 5)
  expect_true(conc3$No_Change == 5)
  expect_true(conc1$Hom_Ref_to_Hom_Alt == 1)
  expect_true(conc2$Hom_Ref_to_Hom_Alt == 0)
  expect_true(conc3$Hom_Ref_to_Hom_Alt == 0)
})

test_that("conta sample swap pairwise genotype concordance", {
  options(warn = -1)
  # Test conta genotype files exists
  expect_true(file.exists(swap_sim_file_1))
  expect_true(file.exists(swap_sim_file_2))
  expect_true(file.exists(swap_sim_file_3))
  # Construct file and label lists
  # files <- list(c(swap_sim_file_1, swap_sim_file_2, swap_sim_file_3))
  # labels <- list(c("sample1", "sample2", "sample3"))

  files <- list(c(swap_sim_file_1, swap_sim_file_2))
  labels <- list(c("sample1", "sample2"))

  # Calculate pairwise concordance, three samples, threshold of 0.7
  concordances <- conta::conta_swap(files, labels, 0.7)
  # Validate the dimensions of the dataframe
  concordances_long <- data.frame(concordances[2])
  expect_true(nrow(concordances_long) == 8)
  expect_true(ncol(concordances_long) == 4)
  # Validate the dimensions of the dataframe
  concordances_wide <- data.frame(concordances[1])
  expect_true(nrow(concordances_wide) == 1)
  expect_true(ncol(concordances_wide) == 10)
  # Validate the expected concordance values
  expect_true(concordances_long %>% filter(
    Sample1 == "sample1" & Sample2 == "sample2" & metric_name == "Concordance") %>%
      select(metric_value) == 0.9)
  # Validate the expected swap calls
  expect_true(concordances_long %>% filter(
    Sample1 == "sample1" & Sample2 == "sample2" & metric_name == "Call") %>%
      select(metric_value) == FALSE)

  # Calculate pairwise concordance, three samples, threshold of 0.4
  concordances <- conta::conta_swap(files, labels, 0.4)
  concordances_long <- data.frame(concordances[2])
  # Validate the dimensions of the dataframe
  expect_true(nrow(concordances_long) == 8)
  expect_true(ncol(concordances_long) == 4)
  # Validate the expected concordance values
  expect_true(concordances_long %>% filter(
    Sample1 == "sample1" & Sample2 == "sample2" & metric_name == "Concordance") %>%
      select(metric_value) == 0.9)
  # Validate the expected swap calls
  expect_true(concordances_long %>% filter(
    Sample1 == "sample1" & Sample2 == "sample2" & metric_name == "Call") %>%
      select(metric_value) == FALSE)


  # Ensure the dimensions of results are correct. Expect that the results are a list of 2.
  concordances <- conta::conta_swap(files, labels, 0.85)
  expect_equal(length(concordances), 2)
})
