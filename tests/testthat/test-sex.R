# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test conta sex calls")

test_that("Test whether sex calling is correct in some cases", {

  chr_y_male_threshold <- 0.0005

  # Clean male
  chr_y_stats <- data.frame(normalized_count = 0.0005)
  result <- data.frame(cf = 0, conta_call = FALSE)
  sex <- conta::get_sex_call(chr_y_stats, chr_y_male_threshold, result)
  expect_true(sex == "Male")

  # Gray zone
  chr_y_stats <- data.frame(normalized_count = 0.0004)
  result <- data.frame(cf = 0, conta_call = FALSE)
  sex <- conta::get_sex_call(chr_y_stats, chr_y_male_threshold, result)
  expect_true(sex == "No_call")

  # Clean female
  chr_y_stats <- data.frame(normalized_count = 0.00001)
  result <- data.frame(cf = 0, conta_call = FALSE)
  sex <- conta::get_sex_call(chr_y_stats, chr_y_male_threshold, result)
  expect_true(sex == "Female")

  # Highly contaminated female
  chr_y_stats <- data.frame(normalized_count = 0.00002)
  result <- data.frame(cf = 0.2, conta_call = TRUE)
  sex <- conta::get_sex_call(chr_y_stats, chr_y_male_threshold, result)
  expect_true(sex == "Female")

  # Highly contaminated male
  chr_y_stats <- data.frame(normalized_count = 0.0005)
  result <- data.frame(cf = 0.2, conta_call = TRUE)
  sex <- conta::get_sex_call(chr_y_stats, chr_y_male_threshold, result)
  expect_true(sex == "No_call")

  # High male contamination
  chr_y_stats <- data.frame(normalized_count = 0.0004)
  result <- data.frame(cf = 0.2, conta_call = TRUE)
  sex <- conta::get_sex_call(chr_y_stats, chr_y_male_threshold, result)
  expect_true(sex == "No_call")

  # Low male contamination
  chr_y_stats <- data.frame(normalized_count = 0.0002)
  result <- data.frame(cf = 0.2, conta_call = TRUE)
  sex <- conta::get_sex_call(chr_y_stats, chr_y_male_threshold, result)
  expect_true(sex == "Female")
})
