context("test loh")


test_that("test loh in wgs setting", {

  # Simulate WGS regions of 10,000 SNPs at depth 50
  er_min <- 0.00003
  er_max <- 0.0001
  dp <- 50
  n <- 10000
  min_maf <- 0.25
  blackswan <- 0.05
  call_thr <- 0.01

  # Simulate het sample, should not be called
  delta <- 0
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont < call_thr)

  # Simulate 5% delta, should be called loh
  delta <- 0.05
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 10% delta, should be called loh
  delta <- 0.1
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 20% delta, should be called loh
  delta <- 0.2
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 50% delta, should be called loh
  delta <- 0.5
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 50% contamination, should be called contamination
  delta <- 0
  alpha <- 0.5
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 20% contamination, should be called contamination
  delta <- 0
  alpha <- 0.2
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 10% contamination, should be called contamination
  delta <- 0
  alpha <- 0.1
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 1% contamination, should be called contamination
  delta <- 0
  alpha <- 0.01
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 0.1% contamination, should be called contamination
  delta <- 0
  alpha <- 0.001
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 20% LOH and 1% contamination, should be called LOH
  delta <- 0.2
  alpha <- 0.01
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont > call_thr & lr_loh > lr_cont)

  # Simulate 10% LOH and 0.1% contamination, should be called LOH
  delta <- 0.1
  alpha <- 0.001
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont > call_thr & lr_loh > lr_cont)

})

test_that("test loh in targeted setting", {

  # Simulate WGS regions of 10,000 SNPs at depth 50
  er_min <- 0.000003
  er_max <- 0.00001
  dp <- 3000
  n <- 20
  min_maf <- 0.25
  blackswan <- 0.05
  call_thr <- 0.01

  # Simulate het sample, should not be called
  delta <- 0
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont < call_thr)

  # Simulate 5% delta, should be called loh
  delta <- 0.05
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 10% delta, should be called loh
  delta <- 0.1
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 20% delta, should be called loh
  delta <- 0.2
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 50% delta, should be called loh
  delta <- 0.5
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont < call_thr & lr_loh > lr_cont)

  # Simulate 50% contamination, should be called contamination
  delta <- 0
  alpha <- 0.5
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 20% contamination, should be called contamination
  delta <- 0
  alpha <- 0.2
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 10% contamination, should be called contamination
  delta <- 0
  alpha <- 0.1
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 1% contamination, should be called contamination
  delta <- 0
  alpha <- 0.01
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 0.1% contamination, should be called contamination
  delta <- 0
  alpha <- 0.001
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_loh < lr_cont)

  # Simulate 20% LOH and 1% contamination, should be called LOH
  delta <- 0.2
  alpha <- 0.01
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont > call_thr & lr_loh > lr_cont)

  # Simulate 10% LOH and 0.1% contamination, should be called LOH
  delta <- 0.1
  alpha <- 0.001
  dat <- simulate_loh_conta(n, min_maf, dp, dp, er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont > call_thr & lr_loh > lr_cont)

})

test_that("test loh in variable length (RNA) setting", {

  # Simulate WGS regions of 10,000 SNPs at depth 50
  er_min <- 0.00003
  er_max <- 0.0001
  dp_min <- 10
  dp_max <- 500000
  n <- 200
  min_maf <- 0.25
  blackswan <- 0.05
  call_thr <- 0.01

  # Simulate het sample, should not be called
  delta <- 0
  alpha <- 0
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont < call_thr)

  # Simulate 10% contaminated sample, should not be called
  delta <- 0
  alpha <- 0.1
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_cont > lr_loh)

  # Simulate 1% contaminated sample, should be called contaminated
  delta <- 0
  alpha <- 0.01
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh < call_thr & lr_cont > call_thr & lr_cont > lr_loh)

  # Simulate 1% contaminated sample with 1% loh, should be called contaminated
  delta <- 0.01
  alpha <- 0.01
  dat <- simulate_loh_conta(n, min_maf, dp_min, dp_max,
                            er_min, er_max, delta, alpha)
  lr_loh <- avg_lr_loh(dat, delta, blackswan)
  lr_cont <- avg_lr_cont(dat, alpha, blackswan)
  expect_true(lr_loh > call_thr & lr_cont > call_thr & lr_cont > lr_loh)
})