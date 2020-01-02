# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test likelihood")

test_that("test likelihood calculation", {

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

test_that("test vectorized likelihood calculation", {

  expect_equal(round(log_lr(c(0.5, 0.5, 0.5),
                            c(30, 30, 30),
                            c(0.1, 0.1, 0.1),
                            c(0.0001, 0.001, 0.0001),
                            c(1, 1, 3)), 1),
               c(3.2, 1.1, 9.3))
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

#' 0.2 * 0.7 / (0.2 * 0.7 + 0.8 * 0.3) = 0.37
#' 0.2 * 0.3 / (0.2 * 0.3 + 0.8 * 0.7) = 0.1
test_that("test allele frequency normalization", {
  expect_true(round(nloh(0.2, 0.2), 2) == 0.37)
  expect_true(round(nloh(0.2, -0.2), 2) == 0.1)
})

# Likelihoods below were verified by alternative implementations of beta
# binomial likelihood in Go
test_that("test beta binomial likelihood", {
  expect_true(round(ldbetabinom(10, 40, .25, 40), 1) == -2.3)
  expect_true(round(dbetabinom(10, 40, .25, 40), 1) == round(exp(-2.3), 1))
  expect_true(round(ldbetabinom(20, 80, .25, 40), 1) == -2.8)
  expect_true(round(dbetabinom(20, 80, .25, 40), 2) == round(exp(-2.8), 2))
  expect_true(round(ldbetabinom(100, 200, .4, 400), 1) == -5.8)
  expect_true(round(dbetabinom(100, 200, .4, 400), 3) == round(exp(-5.8), 3))

})
