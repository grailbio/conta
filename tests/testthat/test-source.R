# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test source detect")

test_that("conta source detection run",  {

  expect_true(file.exists(in_dir_source))
  source_tsvs <- list.files(path = in_dir_source,
                            pattern = "*.nobaq.tsv", full.names = TRUE)
  dir.create(out_dir_source, recursive = TRUE)

  # Run conta on a set titration experiments with known source
  for (tsv_file in source_tsvs) {
    name <- basename(gsub(".subset.collapsed.nobaq.tsv", "", tsv_file))
    maf_tsv <- gsub("nobaq.tsv", "nobaq.maf.tsv", tsv_file)
    maf_tsv <- paste(out_dir_source, name, basename(maf_tsv), sep = "/")
    dir.create(dirname(maf_tsv), recursive = TRUE)
    conta::intersect_snps(tsv_file, maf_tsv, dbSNP_file_art, FALSE)
    conta_main(tsv_file = maf_tsv, sample = name, cores = 8,
               save_dir = dirname(maf_tsv), min_cf = 0.00005, lr_th = 0.05)
  }

  # Run conta source
  conta_source(base = out_dir_source, out_file = out_source, cores = 4)

  expect_true(file.exists(out_source))
  result <- read_data_table(out_source)
  expect_true(result[, .N] == 6)
  expect_false(result[1, source_call])
  expect_equal(sum(result[2:5, source_call]), 4)
  expect_equal(sum(result[2:5, best_sample] == "170410_cfdna_100T_34795B_1"), 4)

  call_mt <- 0.005
  expect_true(result[1, best_gt_lr] < call_mt + result[1, avg_maf_lr])
  expect_true(result[2, best_gt_lr] > call_mt + result[2, avg_maf_lr])
  expect_true(result[3, best_gt_lr] > call_mt + result[3, avg_maf_lr])
  expect_true(result[4, best_gt_lr] > call_mt + result[4, avg_maf_lr])
  expect_true(result[5, best_gt_lr] > call_mt + result[5, avg_maf_lr])
})

test_that("conta source detection run single",  {

  expect_true(file.exists(in_dir_source_single))
  source_tsvs <- list.files(path = in_dir_source_single,
                            pattern = "*.nobaq.tsv", full.names = TRUE)
  dir.create(out_dir_source_single, recursive = TRUE)

  # Run conta on a set titration experiments with known source
  for (tsv_file in source_tsvs) {
    name <- basename(gsub(".subset.collapsed.nobaq.tsv", "", tsv_file))
    maf_tsv <- gsub("nobaq.tsv", "nobaq.maf.tsv", tsv_file)
    maf_tsv <- paste(out_dir_source_single, name, basename(maf_tsv), sep = "/")
    dir.create(dirname(maf_tsv), recursive = TRUE)
    conta::intersect_snps(tsv_file, maf_tsv, dbSNP_file_art, FALSE)
    conta_main(tsv_file = maf_tsv, sample = name, cores = 8,
               save_dir = dirname(maf_tsv), min_cf = 0.00005, lr_th = 0.05)
  }

  # Run conta source
  conta_source(base = out_dir_source_single, out_file = out_source_single, cores = 4)

  expect_true(file.exists(out_source_single))
  result <- read_data_table(out_source_single)
  expect_true(result[, .N] == 1)
  expect_false(result[1, source_call])
  expect_true(is.na(result[1, best_sample]))
  expect_true(is.na(result[1, best_gt_lr]))
})

test_that("conta source detection run double",  {

  expect_true(file.exists(in_dir_source_double))
  source_tsvs <- list.files(path = in_dir_source_double,
                            pattern = "*.nobaq.tsv", full.names = TRUE)
  dir.create(out_dir_source_double, recursive = TRUE)

  # Run conta on a set titration experiments with known source
  for (tsv_file in source_tsvs) {
    name <- basename(gsub(".subset.collapsed.nobaq.tsv", "", tsv_file))
    maf_tsv <- gsub("nobaq.tsv", "nobaq.maf.tsv", tsv_file)
    maf_tsv <- paste(out_dir_source_double, name, basename(maf_tsv), sep = "/")
    dir.create(dirname(maf_tsv), recursive = TRUE)
    conta::intersect_snps(tsv_file, maf_tsv, dbSNP_file_art, FALSE)
    conta_main(tsv_file = maf_tsv, sample = name, cores = 8,
               save_dir = dirname(maf_tsv), min_cf = 0.00005,
               lr_th = 0.05)
  }

  # Also add an empty conta result
  empty_out <- paste(out_dir_source_double, "sample1",
                     "sample1.conta.tsv", sep = "/")
  dir.create(dirname(empty_out))
  utils::write.table(empty_result("sample1"), file = empty_out, sep = "\t",
                     row.names = FALSE, quote = FALSE)

  # Run conta source
  conta_source(base = out_dir_source_double, out_file = out_source_double, cores = 1)

  expect_true(file.exists(out_source_double))
  result <- read_data_table(out_source_double)
  expect_true(result[, .N] == 2)
  expect_false(result[1, source_call])
  expect_true(result[2, source_call])
  expect_true(result[2, best_sample] == "170410_cfdna_100T_34795B_1")
  call_mt <- 0.005
  expect_true(result[1, best_gt_lr] < call_mt + result[1, avg_maf_lr])
  expect_true(result[2, best_gt_lr] > call_mt + result[2, avg_maf_lr])
})
