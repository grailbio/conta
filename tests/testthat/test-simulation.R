context("test simulation")

test_that("conta test run with no simulation on clean sample", {
  
  expect_true(file.exists(sim_file))
  
  # run conta on dummy sample
  conta_main(sim_file, "sim_0perc", out_dir_sim, sim_level = 0, cores = 8)
  conta_out <- file.path(out_dir_sim, "sim_0perc.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  expect_false(result[, conta_call])
})

test_that("conta test run with simulation 0.2 percent", {
  
  expect_true(file.exists(sim_file))
  
  # run conta on dummy sample
  conta_main(sim_file, "sim_0.2perc", out_dir_sim, sim_level = 0.002, cores = 8)
  conta_out <- file.path(out_dir_sim, "sim_0.2perc.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  expect_true(result[, conta_call])
  expect_equal(result[, cf], 0.002, tolerance = 1e-3)
})

test_that("conta test run with simulation 0.5 percent", {
  
  expect_true(file.exists(sim_file))
  
  # run conta on dummy sample
  conta_main(sim_file, "sim_0.5perc", out_dir_sim, sim_level = 0.005, cores = 8)
  conta_out <- file.path(out_dir_sim, "sim_0.5perc.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  expect_true(result[, conta_call])
  expect_equal(result[, cf], 0.005, tolerance = 2e-3)
})

test_that("conta test run with simulation 1 percent", {
  
  expect_true(file.exists(sim_file))
  
  # run conta on dummy sample
  conta_main(sim_file, "sim_1perc", out_dir_sim, sim_level = 0.01, cores = 8)
  conta_out <- file.path(out_dir_sim, "sim_1perc.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  expect_true(result[, conta_call])
  expect_equal(result[, cf], 0.01, tolerance = 2e-3)
})

test_that("conta test run with simulation 5 percent", {
  
  expect_true(file.exists(sim_file))
  
  # run conta on dummy sample
  conta_main(sim_file, "sim_5perc", out_dir_sim, sim_level = 0.05, cores = 8)
  conta_out <- file.path(out_dir_sim, "sim_5perc.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  expect_true(result[, conta_call])
  expect_equal(result[, cf], 0.05, tolerance = 1e-2)
})

test_that("conta test run with simulation 10 percent", {
  
  expect_true(file.exists(sim_file))
  
  # run conta on dummy sample
  conta_main(sim_file, "sim_10perc", out_dir_sim, sim_level = 0.1, cores = 8)
  conta_out <- file.path(out_dir_sim, "sim_10perc.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  expect_true(result[, conta_call])
  expect_equal(result[, cf], 0.1, tolerance = 1e-2)
})

test_that("conta test run with simulation 20 percent", {
  
  expect_true(file.exists(sim_file))
  
  # run conta on dummy sample
  conta_main(sim_file, "sim_20perc", out_dir_sim, sim_level = 0.2, cores = 8)
  conta_out <- file.path(out_dir_sim, "sim_20perc.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  expect_true(result[, conta_call])
  expect_equal(result[, cf], 0.2, tolerance = 3e-2)
})
