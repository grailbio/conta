# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test likelihood")

test_that("test exact likelihood calculation", {

  expect_true(round(log_lr(0.5, 30, 0.1, 0.0001, 1), 1) == 3.2)
  expect_true(round(log_lr(0.5, 30, 0.1, 0.001, 1), 1) == 1.1)
  expect_true(round(log_lr(0.5, 30, 0.1, 0.0001, 3), 1) == 9.3)
  expect_true(round(log_lr(0.5, 30, 0.1, 0.001, 3), 1) == 9)
  expect_true(round(log_lr(0.5, 30, 0.1, 0.0001, 5), 1) == 8.4)
  expect_true(round(log_lr(0.5, 30, 0.1, 0.001, 5), 1) == 8.5)
  expect_true(round(log_lr(0.5, 30, 0.1, 0.0001, 10), 1) == 2.9)
  expect_true(round(log_lr(0.5, 30, 0.1, 0.001, 10), 1) == 2.9)
  expect_true(round(log_lr(0.5, 30, 0.01, 0.0001, 10), 1) == 0)
  expect_true(round(log_lr(0.5, 30, 0.01, 0.001, 10), 1) == 0)
  expect_true(round(log_lr(0.5, 30, 0.01, 0.0001, 5), 1) == 0.4)
  expect_true(round(log_lr(0.5, 30, 0.01, 0.001, 5), 1) == 0.6)
  expect_true(round(log_lr(0.5, 30, 0.01, 0.0001, 3), 1) == 5)
  expect_true(round(log_lr(0.5, 30, 0.01, 0.001, 3), 1) == 4.9)
  expect_true(round(log_lr(0.5, 30, 0.01, 0.0001, 1), 1) == 3.6)
  expect_true(round(log_lr(0.5, 30, 0.01, 0.001, 1), 1) == 1.5)
})


test_that("test likelihood ratio of contamination is higher when AD > 0", {
  expect_true(log_lr(0.5, 30, 0.01, 0.0003, 1) >
              log_lr(0.5, 30, 0.01, 0.0003, 0))
  expect_true(log_lr(0.5, 30, 0.01, 0.0003, 2) >
              log_lr(0.5, 30, 0.01, 0.0003, 0))
  expect_true(log_lr(0.5, 30, 0.01, 0.0003, 3) >
              log_lr(0.5, 30, 0.01, 0.0003, 0))
})

test_that("test likelihood ratio of contamination is higher when
          contamination probability is higher", {
  expect_true(log_lr(0.5, 30, 0.01, 0.0003, 1) >
              log_lr(0.5, 30, 0.01, 0.0003, 0))
  expect_true(log_lr(0.6, 30, 0.01, 0.0003, 1) >
              log_lr(0.5, 30, 0.01, 0.0003, 0))
  expect_true(log_lr(0.7, 30, 0.01, 0.0003, 1) >
              log_lr(0.5, 30, 0.01, 0.0003, 0))
})
