context("test sample swap")

test_that("conta sample swap detection run", {

  options(warn = -1)

  # Test conta genotype files exists
  expect_true(file.exists(swap_test_file_1))
  expect_true(file.exists(swap_test_file_2))
  expect_true(file.exists(swap_test_file_3))

  # Load the genotype files
  dt1 <- load_conta_file(swap_test_file_1, targeted = FALSE)
  dt2 <- load_conta_file(swap_test_file_2, targeted = FALSE)
  dt3 <- load_conta_file(swap_test_file_3, targeted = FALSE)

  # Test the data tables contain expected data
  expect_equal(nrow(dt1), 9999)
  expect_equal(nrow(dt2), 9999)
  expect_equal(nrow(dt3), 9999)

  # Calculate concordance between c1 & c2 and c1 & c3
  conc1 <- genotype_concordance(dt1, dt2)
  conc2 <- genotype_concordance(dt1, dt3)
  conc3 <- genotype_concordance(dt2, dt3)

  expect_true(conc1 > 0.95)
  expect_true(conc2 < 0.7)
  expect_true(conc3 < 0.7)

})
