# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test parsing")

test_that("conta test run with no simulation on clean sample", {
  samples <- paste0("sample", 1:10, "_", rep(c(0, 0.1, 0.5), each = 10))
  avg_log_lr <- rnorm(30)
  test_data <- data.frame(sample = samples,
                          avg_log_lr = avg_log_lr)
  write.table(test_data,
              "test_sim.tsv",
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE)
  got <- read_sim_results("test_sim.tsv")
  expect_equal(dim(got), c(30, 3))
  expect_equal(colnames(got), c("sample", "avg_log_lr", "conta_level"))
  expect_equal(got$sample, rep(paste0("sample", 1:10), 3))
  expect_equal(got$avg_log_lr, avg_log_lr)
  expect_equal(got$conta_level, rep(c(0, 0.1, 0.5), each = 10))

  # expect error if samples are misformatted
  test_data$sample[10] <- NA
  write.table(test_data,
              "test_sim.tsv",
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE)
  expect_error(read_sim_results("test_sim.tsv"))

  test_data$sample[10] <- "sample3_contaLevel_high"
  write.table(test_data,
              "test_sim.tsv",
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE)
  expect_error(read_sim_results("test_sim.tsv"))

  # expect error if missing required columns
  test_data$sample <- NULL
  write.table(test_data,
              "test_sim.tsv",
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE)
  expect_error(read_sim_results("test_sim.tsv"))
  file.remove("test_sim.tsv")
})
