context("test source detect")

test_that("conta source detection run", {
  
  expect_true(file.exists(out_dir))
  source_out <- paste(out_dir, "source.tsv", sep = "/")
  conta_source(out_dir, source_out, cores = 4)
  expect_true(file.exists(source_out))
  result <- read_data_table(source_out)
  expect_true(result[, .N] == 2)
  expect_false(result[1, source_call])
  expect_false(result[2, source_call])
  expect_true(result[1, is.na(best_sample)])
  expect_true(result[2, is.na(best_sample)])
  expect_equal(result[1, best_gt_score], 0)
  expect_equal(result[2, best_gt_score], 0)
})
