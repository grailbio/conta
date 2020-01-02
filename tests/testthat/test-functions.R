# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test functions")

test_that("Test skewed sigmoid", {
  expect_true(ssig(0, 0.5) == 0.5)
  expect_true(ssig(0, 0.3) == 0.3)
  expect_true(ssig(0, 0.7) == 0.7)
  expect_true(ssig(5) == 1 / (1 + exp(-5)))
  expect_true(ssig(-5) == 1 / (1 + exp(5)))
  expect_true(ssig(5, 0.3) < ssig(5))
  expect_true(ssig(5) < ssig(5, 0.7))
  expect_true(ssig(-5, 0.3) < ssig(-5))
  expect_true(ssig(-5) < ssig(-5, 0.7))
})

