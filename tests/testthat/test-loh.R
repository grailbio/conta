# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test LOH")


test_that("test LOH in wgs setting", {

  # Simulate WGS regions of 10,000 SNPs at depth 50
  er_min <- 0.00003
  er_max <- 0.0001
  dp <- 50
  n <- 10000
  min_maf <- 0.25
  blackswan <- 1e-6
  call_thr <- 0.001

  # Simulate het sample, should not be called
  delta <- 0
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont < call_thr)

  # Simulate 5% delta, should be called LOH
  delta <- 0.05
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 10% delta, should be called LOH
  delta <- 0.1
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 20% delta, should be called LOH
  delta <- 0.2
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 50% delta, should be called LOH
  delta <- 0.5
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 50% contamination, should be called contamination
  delta <- 0
  alpha <- 0.5
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 20% contamination, should be called contamination
  delta <- 0
  alpha <- 0.2
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 10% contamination, should be called contamination
  delta <- 0
  alpha <- 0.1
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 1% contamination, should be called contamination
  delta <- 0
  alpha <- 0.01
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 0.1% contamination, should be called contamination
  delta <- 0
  alpha <- 0.001
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 20% LOH and 1% contamination, should be called LOH
  delta <- 0.2
  alpha <- 0.001
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont > call_thr & lr_loh > lr_cont)

  # Simulate 10% LOH and 0.1% contamination, should be called LOH
  delta <- 0.1
  alpha <- 0.001
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont > call_thr & lr_loh > lr_cont)

})

test_that("test LOH in targeted setting", {

  # Simulate targeted regions of 20 SNPs at depth 1000
  er_min <- 0.000003
  er_max <- 0.00001
  dp <- 1000
  n <- 25
  min_maf <- 0.25
  blackswan <- 1e-6
  call_thr <- 0.001

  # Simulate het sample, should not be called
  delta <- 0
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont < call_thr)

  # Simulate 5% delta, should be called LOH
  delta <- 0.05
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 10% delta, should be called LOH
  delta <- 0.1
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 20% delta, should be called LOH
  delta <- 0.2
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 50% delta, should be called LOH
  delta <- 0.5
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 50% contamination, should be called contamination
  delta <- 0
  alpha <- 0.5
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 20% contamination, should be called contamination
  delta <- 0
  alpha <- 0.2
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 10% contamination, should be called contamination
  delta <- 0
  alpha <- 0.1
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 1% contamination, should be called contamination
  delta <- 0
  alpha <- 0.01
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 0.1% contamination, should be called contamination
  delta <- 0
  alpha <- 0.001
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 20% LOH and 1% contamination, should be called LOH
  delta <- 0.2
  alpha <- 0.001
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont > call_thr & lr_loh > lr_cont)

  # Simulate 10% LOH and 0.1% contamination, should be called LOH
  delta <- 0.1
  alpha <- 0.001
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont > call_thr & lr_loh > lr_cont)

})

test_that("test LOH in variable length (RNA) setting", {

  # Simulate RNA regions of 200 SNPs at depth 10 to 10000
  er_min <- 0.00003
  er_max <- 0.0001
  dp_min <- 10
  dp_max <- 10000
  n <- 200
  min_maf <- 0.25
  blackswan <- 1e-6
  call_thr <- 0.005

  # Simulate het sample, should not be called
  delta <- 0
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont < call_thr)

  # Simulate 10% contaminated sample, should not be called
  delta <- 0
  alpha <- 0.1
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_cont > lr_loh)

  # Simulate 20% contaminated sample, should not be called
  delta <- 0
  alpha <- 0.2
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_cont > lr_loh)

  # Simulate 50% contaminated sample, should not be called
  delta <- 0
  alpha <- 0.5
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_cont > lr_loh)

  # Simulate 10% LOH, should be called LOH
  delta <- 0.1
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_cont < lr_loh)

  # Simulate 20% LOH, should be called LOH
  delta <- 0.2
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_cont < lr_loh)

  # Simulate 50% LOH, should be called LOH
  delta <- 0.5
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_cont < lr_loh)


  # Simulate 0.1% contaminated sample, should be called contaminated
  delta <- 0
  alpha <- 0.001
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_cont > lr_loh)

  # Simulate 1% contaminated sample with 1% LOH, should be called contaminated
  delta <- 0.01
  alpha <- 0.01
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_llr_loh(dat, delta, blackswan)
  lr_cont <- avg_llr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont > call_thr & lr_cont > lr_loh)
})
