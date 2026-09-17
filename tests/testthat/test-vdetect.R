library(testthat)
library(enviPat)

data(isotopes)

# Shared fixtures
pat_caff <- isopattern(isotopes, "C8H10N4O2",
                       threshold = 0.001, plotit = FALSE,
                       charge = FALSE, algo = 1, verbose = FALSE)
env_caff <- envelope(pat_caff, ppm = FALSE, dmz = 0.0001,
                     env = "Gaussian", resolution = 1e5,
                     plotit = FALSE, verbose = FALSE)

# ── Output structure ──────────────────────────────────────────────────────────

test_that("vdetect returns named list of matrices with m/z and abundance columns", {
  result <- vdetect(env_caff, detect = "centroid",
                    plotit = FALSE, verbose = FALSE)
  expect_type(result, "list")
  expect_equal(length(result), 1L)
  expect_named(result, "C8H10N4O2")
  expect_true(is.matrix(result[[1]]))
  expect_equal(ncol(result[[1]]), 2L)
  expect_equal(colnames(result[[1]]), c("m/z", "abundance"))
})

# ── detect="centroid" (type=0) ────────────────────────────────────────────────

test_that("centroid m/z values are positive and sorted ascending", {
  result <- vdetect(env_caff, detect = "centroid",
                    plotit = FALSE, verbose = FALSE)
  mz <- result[[1]][, "m/z"]
  expect_true(all(mz > 0))
  expect_true(all(diff(mz) > 0))
})

test_that("centroid max abundance equals 100", {
  result <- vdetect(env_caff, detect = "centroid",
                    plotit = FALSE, verbose = FALSE)
  expect_equal(max(result[[1]][, "abundance"]), 100)
})

test_that("centroid m/z values lie within envelope m/z range", {
  result <- vdetect(env_caff, detect = "centroid",
                    plotit = FALSE, verbose = FALSE)
  env_range <- range(env_caff[[1]][, "m/z"])
  expect_true(all(result[[1]][, "m/z"] >= env_range[1]))
  expect_true(all(result[[1]][, "m/z"] <= env_range[2]))
})

test_that("centroid detects at least one peak", {
  result <- vdetect(env_caff, detect = "centroid",
                    plotit = FALSE, verbose = FALSE)
  expect_gt(nrow(result[[1]]), 0L)
})

# ── detect="intensoid" (type=1) ───────────────────────────────────────────────

test_that("intensoid detects same number of peaks as centroid", {
  cent  <- vdetect(env_caff, detect = "centroid",
                   plotit = FALSE, verbose = FALSE)
  inten <- vdetect(env_caff, detect = "intensoid",
                   plotit = FALSE, verbose = FALSE)
  expect_equal(nrow(cent[[1]]), nrow(inten[[1]]))
})

test_that("intensoid m/z values are positive and sorted ascending", {
  result <- vdetect(env_caff, detect = "intensoid",
                    plotit = FALSE, verbose = FALSE)
  expect_true(all(result[[1]][, "m/z"] > 0))
  expect_true(all(diff(result[[1]][, "m/z"]) > 0))
})

test_that("intensoid max abundance equals 100", {
  result <- vdetect(env_caff, detect = "intensoid",
                    plotit = FALSE, verbose = FALSE)
  expect_equal(max(result[[1]][, "abundance"]), 100)
})

test_that("centroid and intensoid m/z values are not identical", {
  cent  <- vdetect(env_caff, detect = "centroid",
                   plotit = FALSE, verbose = FALSE)
  inten <- vdetect(env_caff, detect = "intensoid",
                   plotit = FALSE, verbose = FALSE)
  expect_false(identical(cent[[1]][, "m/z"], inten[[1]][, "m/z"]))
})

# ── detect="valley" (type=2) ──────────────────────────────────────────────────

test_that("valley m/z values are positive and sorted ascending", {
  result <- vdetect(env_caff, detect = "valley",
                    plotit = FALSE, verbose = FALSE)
  expect_true(all(result[[1]][, "m/z"] > 0))
  expect_true(all(diff(result[[1]][, "m/z"]) > 0))
})

test_that("valley count is at least centroid count minus 1", {
  cent   <- vdetect(env_caff, detect = "centroid",
                    plotit = FALSE, verbose = FALSE)
  valley <- vdetect(env_caff, detect = "valley",
                    plotit = FALSE, verbose = FALSE)
  expect_gte(nrow(valley[[1]]), nrow(cent[[1]]) - 1L)
})

# ── Pipeline cross-check ──────────────────────────────────────────────────────

test_that("centroid m/z range is close to isopattern m/z range", {
  cent      <- vdetect(env_caff, detect = "centroid",
                       plotit = FALSE, verbose = FALSE)
  pat_range  <- range(pat_caff[[1]][, "m/z"])
  cent_range <- range(cent[[1]][, "m/z"])
  expect_gte(cent_range[1], pat_range[1] - 0.1)
  expect_lte(cent_range[2], pat_range[2] + 0.1)
})

# ── Error cases ───────────────────────────────────────────────────────────────

test_that("vdetect stops for invalid detect argument", {
  expect_error(
    vdetect(env_caff, detect = "badtype", plotit = FALSE, verbose = FALSE)
  )
})
