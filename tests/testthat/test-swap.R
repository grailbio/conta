# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' @importFrom data.table fread
#' @importFrom dplyr pull
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
  concordances <- conta::conta_swap(files, labels, 0.7, 0, 0)
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
  concordances <- conta::conta_swap(files, labels, 0.4, 0, 0)
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
  concordances <- conta::conta_swap(files, labels, 0.85, 0, 0)
  expect_equal(length(concordances), 2)
})

test_that("conta gt filtering for swap", {
  # Test that files exist
  expect_true(file.exists(filter_tm_gt))
  expect_true(file.exists(filter_wgbs_gt))

  # Check filtering of methyl_3.16 R conta output
  wgbs_dt <- data.table::fread(filter_wgbs_gt)
  expect_equal(nrow(wgbs_dt), 5, info = wgbs_dt)
  wgbs_filtered <- conta::filter_gt_file(wgbs_dt, 0.01, 0)
  expect_equal(nrow(wgbs_filtered), 3, info = wgbs_filtered)
  wgbs_filtered <- conta::filter_gt_file(wgbs_dt, 0.01, 15)
  expect_equal(nrow(wgbs_filtered), 2, info = wgbs_filtered)
  wgbs_filtered <- conta::filter_gt_file(wgbs_dt, 0, 0)
  expect_equal(nrow(wgbs_filtered), 4, info = wgbs_filtered)
  expect_equal(names(wgbs_filtered), c("rsid", "cp", "dp", "er", "gt", "vr", "maf"))
  expect_true(all(startsWith(pull(wgbs_filtered, rsid), "rs")), info = wgbs_filtered)
  expect_true(all(!startsWith(pull(wgbs_filtered, rsid), "rsrs")), info = wgbs_filtered)

  # Check filtering of methyl_3.16 go conta output
  tm_dt <- data.table::fread(filter_tm_gt)
  expect_equal(nrow(tm_dt), 4, info = tm_dt)
  tm_filtered <- conta::filter_gt_file(tm_dt, 0.01, 0)
  expect_equal(nrow(tm_filtered), 1, info = tm_filtered)
  tm_filtered <- conta::filter_gt_file(tm_dt, 0.01, 15)
  expect_equal(nrow(tm_filtered), 0, info = tm_filtered)
  tm_filtered <- conta::filter_gt_file(tm_dt, 0, 0)
  expect_equal(nrow(tm_filtered), 2, info = tm_filtered)
  expect_equal(names(tm_filtered), c("rsid", "cp", "dp", "er", "gt", "vr", "maf", "loh"))
  expect_true(all(startsWith(pull(tm_filtered, rsid), "rs")), info = tm_filtered)
  expect_true(all(!startsWith(pull(tm_filtered, rsid), "rsrs")), info = tm_filtered)
  expect_true(all(!pull(tm_filtered, loh)), info = tm_filtered)

  # Check that concordance can be calculated between them
  concordance <- conta::genotype_concordance(wgbs_filtered, tm_filtered)
  expect_equal(pull(concordance, Concordance), 1, info = concordance)
})
