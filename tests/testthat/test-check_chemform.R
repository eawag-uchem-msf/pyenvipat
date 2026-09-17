library(testthat)
library(enviPat)

data(isotopes)

test_that("check_chemform returns data.frame with correct columns", {
  result <- check_chemform(isotopes, "C2H6O1")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("warning", "new_formula", "monoisotopic_mass"))
  expect_equal(nrow(result), 1L)
})

test_that("check_chemform valid formulas have warning=FALSE", {
  result <- check_chemform(isotopes, c("C2H6O1", "C8H10N4O2"))
  expect_false(result[1, "warning"])
  expect_false(result[2, "warning"])
})

test_that("check_chemform returns correct monoisotopic mass for ethanol", {
  result <- check_chemform(isotopes, "C2H6O1")
  # C2H6O: 2*12 + 6*1.0078250319 + 15.9949146221 = 46.0418...
  expect_equal(result[1, "monoisotopic_mass"], 46.04186, tolerance = 1e-4)
})

test_that("check_chemform returns correct monoisotopic mass for caffeine", {
  result <- check_chemform(isotopes, "C8H10N4O2")
  expect_equal(result[1, "monoisotopic_mass"], 194.0804, tolerance = 1e-4)
})

test_that("check_chemform invalid formula gives warning=TRUE and mass=-9999", {
  result <- check_chemform(isotopes, "ZZZZ999")
  expect_true(result[1, "warning"])
  expect_equal(result[1, "monoisotopic_mass"], -9999)
})

test_that("check_chemform handles mixed valid and invalid formulas", {
  result <- check_chemform(isotopes, c("C2H6O1", "ZZZZ999"))
  expect_false(result[1, "warning"])
  expect_true(result[2, "warning"])
})

test_that("check_chemform get_list=TRUE returns named list with element counts", {
  result <- check_chemform(isotopes, "C8H10N4O2", get_list = TRUE)
  expect_type(result, "list")
  expect_equal(length(result), 1L)
  elem <- result[[1]]
  expect_equal(elem[["C"]], 8)
  expect_equal(elem[["H"]], 10)
  expect_equal(elem[["N"]], 4)
  expect_equal(elem[["O"]], 2)
})

test_that("check_chemform accepts isotope-labeled formula", {
  result <- check_chemform(isotopes, "[13]C2C6H10N4O2")
  expect_false(result[1, "warning"])
  # [13]C adds ~2 Da (2 heavy carbons) → mass > 194
  expect_gt(result[1, "monoisotopic_mass"], 194)
})
