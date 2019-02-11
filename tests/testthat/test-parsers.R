# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test parsing")

test_that("SNP blacklist is parsed and applied", {

  expect_true(file.exists(error_tsv))
  dat1 <- read_and_prep(error_tsv)

  # Make a dummy baseline
  baseline <- data.frame(rsid = sample(dat1$rsid, 5),
                         blacklist = sample(c(TRUE, FALSE),
                                            5, replace = TRUE))

  dat2 <- read_and_prep(error_tsv, baseline = baseline)

  expect_equal(nrow(dat1) - sum(baseline$blacklist),
               nrow(dat2))
})

test_that("SNP blacklist is malformed", {

  expect_true(file.exists(error_tsv))
  dat1 <- read_and_prep(error_tsv)

  # Make a dummy baseline
  baseline <- data.frame(id = sample(dat1$rsid, 5),
                         blaklis = sample(c(TRUE, FALSE),
                                            5, replace = TRUE))

  expect_warning(dat2 <- read_and_prep(error_tsv, baseline = baseline))

  expect_equal(nrow(dat1), nrow(dat2))
})

test_that("Baseline is read from file and blacklist applied", {
  expect_true(file.exists(error_tsv))
  expect_true(file.exists(baseline_tsv))
  dat1 <- read_and_prep(error_tsv)
  baseline <- read_data_table(baseline_tsv)
  dat2 <- read_and_prep(error_tsv, baseline = baseline)
  expect_equal(nrow(dat1) - sum(baseline$blacklist),
               nrow(dat2))
})

test_that("Baseline file specified does not exist.", {
  expect_true(file.exists(error_tsv))
  dat1 <- read_and_prep(error_tsv)
  expect_warning(baseline <- read_data_table(paste(baseline_tsv, ".",
                                                   sep = "")))
  dat2 <- read_and_prep(error_tsv, baseline = baseline)
  expect_equal(nrow(dat1), nrow(dat2))
})

test_that("conta test run with no simulation on clean sample", {
  samples <- paste0("sample", 1:10, "_", rep(c(0, 0.1, 0.5), each = 10))
  avg_log_lr <- rnorm(30)
  test_data <- data.frame(sample = samples,
                          avg_log_lr = avg_log_lr,
                          stringsAsFactors = FALSE)
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

  samples <- paste0("this_is_a_test", 1:10, "_", rep(c(0, 0.1, 0.5), each = 10))
  avg_log_lr <- rnorm(30)
  test_data <- data.frame(sample = samples,
                          avg_log_lr = avg_log_lr,
                          stringsAsFactors = FALSE)
  write.table(test_data,
              "test_sim.tsv",
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE)
  got <- read_sim_results("test_sim.tsv")
  expect_equal(dim(got), c(30, 3))
  expect_equal(colnames(got), c("sample", "avg_log_lr", "conta_level"))
  expect_equal(got$sample, rep(paste0("this_is_a_test", 1:10), 3))
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

test_that("bayesian genotype concordance with hard-cutoffs", {
  dat <- readr::read_tsv(dat_tsv)
  dat <- conta:::bayesian_genotype(dat)
  expect_true(sum(dat$bayes_gt == dat$gt)/nrow(dat) > 0.98)
})

test_that("specify precision testing", {
  y <- c("1.1234", "-1.1234", "13123.1111")
  expected <- c("1.123", "-1.123", "13123.111")
  actual <- conta:::specify_precision(y, 3)
  expect_true(all.equal(actual, expected))
})

test_that("bayesian genotype boundary testing", {

  # Testing upper minor ratio boundary with multiple error rates
  test1 <- data.frame(major_count = 10, minor_count = 10:200, er = 0.0002)
  test1 <- test1 %>% mutate(minor_ratio = minor_count/(major_count+minor_count))
  test1_out <- conta:::bayesian_genotype(test1)
  # Finding the minimum minor_ratio at which the genotype is called as hom_alt
  testthat::expect_true(all.equal(test1_out[
    min(which(test1_out$bayes_gt == "1/1")), ]$minor_ratio, 0.917, 0.001))

  test2 <- data.frame(major_count = 10, minor_count = 10:200, er = 0.0001)
  test2 <- test2 %>% mutate(minor_ratio = minor_count/(major_count+minor_count))
  test2_out <- conta:::bayesian_genotype(test2)
  # Finding the minimum minor_ratio at which the genotype is called as hom_alt
  testthat::expect_true(all.equal(test2_out[
    min(which(test2_out$bayes_gt == "1/1")), ]$minor_ratio, 0.924, 0.001))

  test3 <- data.frame(major_count = 10, minor_count = 10:200, er = 0.00001)
  test3 <- test3 %>% mutate(minor_ratio = minor_count/(major_count+minor_count))
  test3_out <- conta:::bayesian_genotype(test3)
  # Finding the minimum minor_ratio at which the genotype is called as hom_alt
  testthat::expect_true(all.equal(test3_out[
    min(which(test3_out$bayes_gt == "1/1")), ]$minor_ratio, 0.939, 0.001))

  # Testing lower minor ratio boundary with multiple error rates
  test4 <- data.frame(major_count = 10:200, minor_count = 10, er = 0.0002)
  test4 <- test4 %>% mutate(minor_ratio = minor_count/(major_count+minor_count))
  test4_out <- conta:::bayesian_genotype(test4)
  # Finding the minimum minor_ratio at which the genotype is called as hom_ref
  testthat::expect_true(all.equal(test4_out[
    min(which(test4_out$bayes_gt == "0/0")), ]$minor_ratio, 0.083, 0.05))

  test5 <- data.frame(major_count = 10:200, minor_count = 10, er = 0.0001)
  test5 <- test5 %>% mutate(minor_ratio = minor_count/(major_count+minor_count))
  test5_out <- conta:::bayesian_genotype(test5)
  # Finding the minimum minor_ratio at which the genotype is called as hom_ref
  testthat::expect_true(all.equal(test5_out[
    min(which(test5_out$bayes_gt == "0/0")), ]$minor_ratio, 0.076, 0.05))

  test6 <- data.frame(major_count = 10:200, minor_count  = 10, er = 0.00001)
  test6 <- test6 %>% mutate(minor_ratio = minor_count/(major_count+minor_count))
  test6_out <- conta:::bayesian_genotype(test6)
  # Finding the minimum minor_ratio at which the genotype is called as hom_ref
  testthat::expect_true(all.equal(test6_out[
    min(which(test6_out$bayes_gt == "0/0")), ]$minor_ratio, 0.061, 0.05))
})
