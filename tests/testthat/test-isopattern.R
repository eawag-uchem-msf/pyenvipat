library(testthat)
library(enviPat)

data(isotopes)

# ── Output structure ──────────────────────────────────────────────────────────

test_that("isopattern returns named list of matrices", {
  result <- isopattern(isotopes, c("C2H6O1", "C8H10N4O2"),
                       threshold = 0.001, plotit = FALSE,
                       charge = FALSE, algo = 1, verbose = FALSE)
  expect_type(result, "list")
  expect_equal(length(result), 2L)
  expect_named(result, c("C2H6O1", "C8H10N4O2"))
  expect_true(is.matrix(result[[1]]))
  expect_true(is.matrix(result[[2]]))
})

test_that("isopattern matrix first two columns are m/z and abundance", {
  result <- isopattern(isotopes, "C8H10N4O2",
                       threshold = 0.001, plotit = FALSE,
                       charge = FALSE, algo = 1, verbose = FALSE)
  expect_equal(colnames(result[[1]])[1], "m/z")
  expect_equal(colnames(result[[1]])[2], "abundance")
})

test_that("isopattern matrix has isotopologue composition columns beyond first two", {
  result <- isopattern(isotopes, "C8H10N4O2",
                       threshold = 0.001, plotit = FALSE,
                       charge = FALSE, algo = 1, verbose = FALSE)
  expect_gt(ncol(result[[1]]), 2L)
  expect_true("12C" %in% colnames(result[[1]]))
})

# ── Physical correctness ──────────────────────────────────────────────────────

test_that("isopattern m/z values are positive and sorted ascending", {
  result <- isopattern(isotopes, "C8H10N4O2",
                       threshold = 0.001, plotit = FALSE,
                       charge = FALSE, algo = 1, verbose = FALSE)
  mz <- result[[1]][, "m/z"]
  expect_true(all(mz > 0))
  expect_true(all(diff(mz) >= 0))
})

test_that("isopattern abundances are in [0, 100] with max = 100", {
  result <- isopattern(isotopes, "C8H10N4O2",
                       threshold = 0.001, plotit = FALSE,
                       charge = FALSE, algo = 1, verbose = FALSE)
  abund <- result[[1]][, "abundance"]
  expect_true(all(abund >= 0))
  expect_true(all(abund <= 100))
  expect_equal(max(abund), 100)
})

test_that("isopattern monoisotopic m/z matches check_chemform mass", {
  checked <- check_chemform(isotopes, "C2H6O1")
  result  <- isopattern(isotopes, "C2H6O1",
                        threshold = 0.001, plotit = FALSE,
                        charge = FALSE, algo = 1, verbose = FALSE)
  expect_equal(unname(result[[1]][1, "m/z"]),
               unname(checked[1, "monoisotopic_mass"]),
               tolerance = 1e-5)
})

test_that("isopattern 12C count in monoisotopic row equals formula carbon count", {
  result <- isopattern(isotopes, "C8H10N4O2",
                       threshold = 0.001, plotit = FALSE,
                       charge = FALSE, algo = 1, verbose = FALSE)
  expect_equal(unname(result[[1]][1, "12C"]), 8)
})

test_that("isopattern caffeine produces more peaks than ethanol", {
  eth <- isopattern(isotopes, "C2H6O1",
                    threshold = 0.001, plotit = FALSE,
                    charge = FALSE, algo = 1, verbose = FALSE)
  caf <- isopattern(isotopes, "C8H10N4O2",
                    threshold = 0.001, plotit = FALSE,
                    charge = FALSE, algo = 1, verbose = FALSE)
  expect_gt(nrow(caf[[1]]), nrow(eth[[1]]))
})

# ── algo=1 vs algo=2 ──────────────────────────────────────────────────────────

test_that("algo=1 and algo=2 produce same row and column count", {
  r1 <- isopattern(isotopes, "C8H10N4O2",
                   threshold = 0.001, plotit = FALSE,
                   charge = FALSE, algo = 1, verbose = FALSE)
  r2 <- isopattern(isotopes, "C8H10N4O2",
                   threshold = 0.001, plotit = FALSE,
                   charge = FALSE, algo = 2, verbose = FALSE)
  expect_equal(nrow(r1[[1]]), nrow(r2[[1]]))
  expect_equal(ncol(r1[[1]]), ncol(r2[[1]]))
})

test_that("algo=1 and algo=2 produce identical m/z and abundances after sorting", {
  r1 <- isopattern(isotopes, "C8H10N4O2",
                   threshold = 0.001, plotit = FALSE,
                   charge = FALSE, algo = 1, verbose = FALSE)
  r2 <- isopattern(isotopes, "C8H10N4O2",
                   threshold = 0.001, plotit = FALSE,
                   charge = FALSE, algo = 2, verbose = FALSE)
  p1 <- r1[[1]][order(r1[[1]][, "m/z"]), ]
  p2 <- r2[[1]][order(r2[[1]][, "m/z"]), ]
  expect_equal(p1[, "m/z"],       p2[, "m/z"],       tolerance = 1e-8)
  expect_equal(p1[, "abundance"], p2[, "abundance"],  tolerance = 1e-6)
})

# ── Charge tests ──────────────────────────────────────────────────────────────

test_that("charge=1 shifts m/z by minus electron mass", {
  emass <- 0.00054857990924
  r_nc  <- isopattern(isotopes, "C8H10N4O2",
                      threshold = 0.001, plotit = FALSE,
                      charge = FALSE, algo = 1, verbose = FALSE)
  r_z1  <- isopattern(isotopes, "C8H10N4O2",
                      threshold = 0.001, plotit = FALSE,
                      charge = 1, algo = 1, verbose = FALSE)
  expect_equal(r_z1[[1]][1, "m/z"],
               r_nc[[1]][1, "m/z"] - emass,
               tolerance = 1e-7)
})

test_that("charge=2 halves m/z minus two electron masses", {
  emass <- 0.00054857990924
  r_nc  <- isopattern(isotopes, "C8H10N4O2",
                      threshold = 0.001, plotit = FALSE,
                      charge = FALSE, algo = 1, verbose = FALSE)
  r_z2  <- isopattern(isotopes, "C8H10N4O2",
                      threshold = 0.001, plotit = FALSE,
                      charge = 2, algo = 1, verbose = FALSE)
  expected <- (r_nc[[1]][1, "m/z"] - 2 * emass) / 2
  expect_equal(r_z2[[1]][1, "m/z"], expected, tolerance = 1e-7)
})

# ── Error cases ───────────────────────────────────────────────────────────────

test_that("isopattern returns 'error' string for unknown formula", {
  result <- isopattern(isotopes, "ZZZZ999",
                       threshold = 0.001, plotit = FALSE,
                       charge = FALSE, algo = 1, verbose = FALSE)
  expect_equal(result[[1]], "error")
})

test_that("isopattern stops for threshold > 100", {
  expect_error(
    isopattern(isotopes, "C2H6O1", threshold = 200, plotit = FALSE, verbose = FALSE)
  )
})

test_that("isopattern stops for invalid algo", {
  expect_error(
    isopattern(isotopes, "C2H6O1", algo = 3, plotit = FALSE, verbose = FALSE)
  )
})
