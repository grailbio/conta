context("test tsv and vcf intersect")

test_that("Test intersection of tsv and dbSnp to generate tsv", {

  conta::intersect(tsvFile, mafFile, dbSnpFile, FALSE)

  expect_true(file.exists(mafFile))

  maf <- read.table(mafFile, header = T, stringsAsFactors = F)
  expect_true(nrow(maf) == 3)
  expect_true(maf$maf[1] == round(0.01078, 4))
  expect_true(maf$maf[2] == round(0.0003994, 4))
  expect_true(maf$maf[3] == round(0.02276, 4))
})

test_that("Test intersection of pileup and dbSnp to generate tsv", {
  
  conta::intersect(pileupFile, mafFile2, dbSnpFile, FALSE)
  
  expect_true(file.exists(mafFile2))
  
  maf <- read.table(mafFile2, header = T, stringsAsFactors = F)
  expect_true(nrow(maf) == 3)
  expect_true(maf$maf[1] == round(0.01078, 4))
  expect_true(maf$maf[2] == round(0.0003994, 4))
  expect_true(maf$maf[3] == round(0.02276, 4))
})

test_that("Test intersection of pileup with dbSnp to supplementary positions", {
  
  conta::intersect(pileupFile, mafFile3, dbSnpFile, TRUE)
  
  expect_true(file.exists(mafFile3))
  
  maf <- read.table(mafFile3, header = T, stringsAsFactors = F)
  expect_true(nrow(maf) == 6)
  expect_true(maf$maf[1] == 0)
  expect_true(maf$maf[2] == round(0.01078, 4))
  expect_true(maf$maf[3] == round(0.0003994, 4))
  expect_true(maf$maf[4] == 0)
  expect_true(maf$maf[5] == round(0.02276, 4))
  expect_true(maf$maf[6] == 0)
})