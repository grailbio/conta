# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test functions for setting thresholds")

test_that("test filtering data", {
  # Make fake data with 12 samples and 3 contamination levels,
  # expect the filtering to only occur based off the avg_log_lr of
  # the 0 conta_level data (i.e first line of avg_log_lr values)
  sim_data  <- data.frame(sample = paste0("sample", 1:12),
                          conta_level = rep(c(0, 0.1, 0.5), each = 12),
                          avg_log_lr = c(rep(0.1, 9), 1.5, 2, 2,
                                         rnorm(12, mean = 2),
                                         rnorm(12, mean = 3)))
  # Expect samples 11 and 12 to be filtered out by the extreme level filter
  # and sample 10 to be filtered by the quantile filter, no samples filtered
  # by MAD filter
  filtered <- filter_data(sim_data)
  expect_equal(filtered, dplyr::filter(sim_data, sample %in% paste0("sample", 1:9)))

  # Now the two samples tied for top are filtered due to
  # quantile filter, sample 10 filtered by MAD filter
  filtered <- filter_data(sim_data, extreme_level = 2.5)
  expect_equal(filtered, dplyr::filter(sim_data, sample %in% paste0("sample", 1:9)))

  # Now expect sample 9 to be filtered by MAD filter
  sim_data$avg_log_lr[9] <- 1.0
  filtered <- filter_data(sim_data)
  expect_equal(filtered, dplyr::filter(sim_data, sample %in% paste0("sample", 1:8)))

  # Expect error if NAs in avg_log_lr
  sim_data$avg_log_lr[10] <- NA
  expect_error(filtered <- filter_data(sim_data, extreme_level = 2.5),
                 regex = "missing values and NaN's not allowed if 'na.rm' is FALSE")

  # Expect warning if no samples with zero simulated contamination
  sim_data$conta_level <- 0.1
  expect_warning(filtered <- filter_data(sim_data, extreme_level = 2.5),
                 regex = "WARNING: Can't filter data without 0 contamination level")
  expect_equal(filtered, sim_data)

})

test_that("test sensitivity and specificity by contamination level", {
  # Make fake data with 12 samples and 3 contamination levels
  sim_data  <- data.frame(sample = paste0("sample", 1:12),
                          conta_level = rep(c(0, 0.1, 0.5), each = 12),
                          avg_log_lr = c(rep(0.1, 8), 1.0, 1.5, 2, 2,
                                         1.0, 1.5, rep(2.5, 10),
                                         rep(2.5, 12)))
  sens_specs <- sens_spec_by_level(sim_data, target_specs = c(0.75, 0.99))
  expect_equal(dim(sens_specs), c(4, 4)) # expect 2 * 2 rows, num conta levels * num specs
  expect_equal(sens_specs$conta_level, c(0.1, 0.1, 0.5, 0.5))
  # expect 1.0 specificity instead of 0.99, as 1.0 is the closest specificity above 0.99
  # that we can achieve with 12 samples
  expect_equal(sens_specs$specificity, c(0.75, 1.0, 0.75, 1.0))
  expect_equal(sens_specs$sensitivity, c(11/12, 10/12, 1.0, 1.0))

  # expect error without 0 contamination level data or with only 0 contamination data
  sim_data$conta_level <- 0.1
  expect_error(sens_specs <- sens_spec_by_level(sim_data),
               regex = "ERROR: Can't get sens/spec without 0 contamination level data")
  sim_data$conta_level <- 0
  expect_error(sens_specs <- sens_spec_by_level(sim_data),
               regex = "ERROR: Can't get sens/spec with only 0 contamination level data")
})
