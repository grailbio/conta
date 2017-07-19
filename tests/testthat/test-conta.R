context("test conta")

test_that("conta test wgs run with cna and bincounts file", {

  expect_true(file.exists(wgs_tsv))
  expect_true(file.exists(bin_file))
  expect_true(file.exists(cna_file))
  
  # run conta on dummy sample
  conta_main(wgs_tsv, "wgs", out_dir_wgs, bin_file, cna_file, cores = 4)
  conta_out <- file.path(out_dir_wgs, "wgs.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  expect_true(result[, conta_call])
  expect_equal(result[, cf], 0.2, tolerance = 5e-2)
  expect_equal(result[, y_count], 0.0021, tolerance = 1e-4)
  expect_equal(result[, final_stein], 0.35, tolerance = 1e-2)
  expect_equal(result[, final_mapd], 0.01, tolerance = 1e-2)
})

test_that("conta test targeted run", {

  expect_true(file.exists(targeted_tsv))

  # run conta on dummy sample
  conta_main(targeted_tsv, "targeted", out_dir_targeted)
  conta_out <- file.path(out_dir_targeted, "targeted.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  expect_true(result[, conta_call])
  expect_equal(result[, cf], 0.001, tolerance = 1e-4)
})
