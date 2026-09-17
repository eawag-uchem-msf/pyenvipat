library(testthat)
library(enviPat)

data(isotopes)

# Shared fixture
pat_caff <- isopattern(isotopes, "C8H10N4O2",
                       threshold = 0.001, plotit = FALSE,
                       charge = FALSE, algo = 1, verbose = FALSE)

# ── Output structure ──────────────────────────────────────────────────────────

test_that("envelope returns named list of matrices with m/z and abundance columns", {
  result <- envelope(pat_caff, ppm = FALSE, dmz = 0.0001,
                     env = "Gaussian", resolution = 1e5,
                     plotit = FALSE, verbose = FALSE)
  expect_type(result, "list")
  expect_equal(length(result), 1L)
  expect_named(result, "C8H10N4O2")
  expect_true(is.matrix(result[[1]]))
  expect_equal(ncol(result[[1]]), 2L)
  expect_equal(colnames(result[[1]]), c("m/z", "abundance"))
})

# ── Physical correctness ──────────────────────────────────────────────────────

test_that("envelope m/z values are strictly ascending", {
  result <- envelope(pat_caff, ppm = FALSE, dmz = 0.0001,
                     env = "Gaussian", resolution = 1e5,
                     plotit = FALSE, verbose = FALSE)
  expect_true(all(diff(result[[1]][, "m/z"]) > 0))
})

test_that("envelope abundances are non-negative with max = 100", {
  result <- envelope(pat_caff, ppm = FALSE, dmz = 0.0001,
                     env = "Gaussian", resolution = 1e5,
                     plotit = FALSE, verbose = FALSE)
  abund <- result[[1]][, "abundance"]
  expect_true(all(abund >= 0))
  expect_equal(max(abund), 100)
})

test_that("envelope m/z range encompasses the isopattern m/z range", {
  result <- envelope(pat_caff, ppm = FALSE, dmz = 0.0001,
                     env = "Gaussian", resolution = 1e5,
                     plotit = FALSE, verbose = FALSE)
  pat_range <- range(pat_caff[[1]][, "m/z"])
  env_range <- range(result[[1]][, "m/z"])
  expect_lte(env_range[1], pat_range[1])
  expect_gte(env_range[2], pat_range[2])
})

# ── ppm=TRUE exercises iso_ppm_Call ──────────────────────────────────────────

test_that("envelope ppm=TRUE produces valid ascending m/z matrix", {
  result <- envelope(pat_caff, ppm = TRUE, dmz = 1.0,
                     env = "Gaussian", resolution = 1e5,
                     plotit = FALSE, verbose = FALSE)
  expect_true(is.matrix(result[[1]]))
  expect_equal(colnames(result[[1]]), c("m/z", "abundance"))
  expect_true(all(result[[1]][, "m/z"] > 0))
  expect_true(all(diff(result[[1]][, "m/z"]) > 0))
})

test_that("envelope ppm=TRUE finer dmz produces more sticks", {
  r_fine   <- envelope(pat_caff, ppm = TRUE, dmz = 0.1,
                       env = "Gaussian", resolution = 1e5,
                       plotit = FALSE, verbose = FALSE)
  r_coarse <- envelope(pat_caff, ppm = TRUE, dmz = 5.0,
                       env = "Gaussian", resolution = 1e5,
                       plotit = FALSE, verbose = FALSE)
  expect_gt(nrow(r_fine[[1]]), nrow(r_coarse[[1]]))
})

# ── Gaussian vs CauchyLorentz ─────────────────────────────────────────────────

test_that("Gaussian and CauchyLorentz both return valid matrices", {
  r_g <- envelope(pat_caff, ppm = FALSE, dmz = 0.0001,
                  env = "Gaussian",      resolution = 1e5,
                  plotit = FALSE, verbose = FALSE)
  r_c <- envelope(pat_caff, ppm = FALSE, dmz = 0.0001,
                  env = "CauchyLorentz", resolution = 1e5,
                  plotit = FALSE, verbose = FALSE)
  expect_true(is.matrix(r_g[[1]]))
  expect_true(is.matrix(r_c[[1]]))
  expect_equal(nrow(r_g[[1]]), nrow(r_c[[1]]))
})

test_that("Gaussian and CauchyLorentz abundances differ", {
  r_g <- envelope(pat_caff, ppm = FALSE, dmz = 0.0001,
                  env = "Gaussian",      resolution = 1e5,
                  plotit = FALSE, verbose = FALSE)
  r_c <- envelope(pat_caff, ppm = FALSE, dmz = 0.0001,
                  env = "CauchyLorentz", resolution = 1e5,
                  plotit = FALSE, verbose = FALSE)
  expect_false(all(r_g[[1]][, "abundance"] == r_c[[1]][, "abundance"]))
})

# ── Batch ─────────────────────────────────────────────────────────────────────

test_that("envelope processes a two-formula batch", {
  pat_batch <- isopattern(isotopes, c("C2H6O1", "C8H10N4O2"),
                          threshold = 0.001, plotit = FALSE,
                          charge = FALSE, algo = 1, verbose = FALSE)
  result <- envelope(pat_batch, ppm = FALSE, dmz = 0.0001,
                     env = "Gaussian", resolution = 1e5,
                     plotit = FALSE, verbose = FALSE)
  expect_equal(length(result), 2L)
  expect_named(result, c("C2H6O1", "C8H10N4O2"))
  expect_true(all(sapply(result, is.matrix)))
})

# ── Error cases ───────────────────────────────────────────────────────────────

test_that("envelope stops for non-list pattern input", {
  expect_error(envelope("not_a_list", plotit = FALSE, verbose = FALSE))
})

test_that("envelope stops for ppm=TRUE with dmz='get'", {
  expect_error(
    envelope(pat_caff, ppm = TRUE, dmz = "get",
             env = "Gaussian", resolution = 1e5,
             plotit = FALSE, verbose = FALSE)
  )
})

test_that("envelope stops for invalid env argument", {
  expect_error(
    envelope(pat_caff, ppm = FALSE, dmz = 0.0001,
             env = "InvalidShape", resolution = 1e5,
             plotit = FALSE, verbose = FALSE)
  )
})
