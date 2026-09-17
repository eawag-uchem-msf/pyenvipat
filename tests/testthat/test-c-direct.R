library(testthat)
library(enviPat)

# Shared fixtures loaded once at module level
data(isotopes)
iso_elem <- as.character(isotopes[, 1])
iso_iso  <- as.character(isotopes[, 2])
iso_mass <- as.numeric(isotopes[, 3])
iso_abu  <- as.numeric(isotopes[, 4])

# Shared profile inputs
.peak_masses <- 194.08
.peak_abund  <- 100.0
.trace       <- seq(193.5, 194.6, by = 0.0001)

# ── iso_ppm_Call ─────────────────────────────────────────────────────────────

test_that("iso_ppm_Call returns a double vector spanning [start, end]", {
  result <- .Call("iso_ppm_Call", 100.0, 101.0, 5.0, PACKAGE = "enviPat")
  expect_type(result, "double")
  expect_gt(length(result), 0L)
  expect_gte(result[1], 100.0)
  # The trace brackets `end`: the loop appends one point past `end` so the
  # range fully covers it. The last point is >= end and overshoots by at
  # most one local ppm step.
  n <- length(result)
  expect_lte(result[n - 1], 101.0)
  expect_gte(result[n], 101.0)
  expect_lte(result[n] - 101.0, result[n] * 5.0 / 1e6)
})

test_that("iso_ppm_Call spacing is proportional to local m/z (true ppm steps)", {
  result <- .Call("iso_ppm_Call", 100.0, 101.0, 5.0, PACKAGE = "enviPat")
  gaps <- diff(result)
  local_ppm <- gaps / result[-length(result)] * 1e6
  expect_true(all(abs(local_ppm - 5.0) < 0.01))
})

test_that("iso_ppm_Call finer ppm produces more points", {
  r_fine   <- .Call("iso_ppm_Call", 100.0, 101.0, 1.0, PACKAGE = "enviPat")
  r_coarse <- .Call("iso_ppm_Call", 100.0, 101.0, 5.0, PACKAGE = "enviPat")
  expect_gt(length(r_fine), length(r_coarse))
})

test_that("iso_ppm_Call wider range produces more points", {
  r_wide   <- .Call("iso_ppm_Call", 100.0, 110.0, 5.0, PACKAGE = "enviPat")
  r_narrow <- .Call("iso_ppm_Call", 100.0, 101.0, 5.0, PACKAGE = "enviPat")
  expect_gt(length(r_wide), length(r_narrow))
})

test_that("iso_ppm_Call degenerate range (start == end) does not crash", {
  r_single <- .Call("iso_ppm_Call", 100.0, 100.0, 5.0, PACKAGE = "enviPat")
  expect_lte(length(r_single), 1L)
})

# ── iso_centroid_Call ─────────────────────────────────────────────────────────

test_that("iso_centroid_Call centroid of synthetic Gaussian is at known centre", {
  mz    <- seq(199.9, 200.1, by = 0.0001)
  abund <- 100 * exp(-0.5 * ((mz - 200.0) / 0.002)^2)
  out   <- .Call("iso_centroid_Call", as.double(mz), as.double(abund),
                 as.integer(0), PACKAGE = "enviPat")
  expect_type(out, "list")
  expect_length(out, 2L)
  expect_length(out[[1]], 1L)
  expect_equal(out[[1]][1], 200.0, tolerance = 1e-3)
  expect_gt(out[[2]][1], 0)
})

test_that("iso_centroid_Call intensoid of synthetic Gaussian is at known centre", {
  mz    <- seq(199.9, 200.1, by = 0.0001)
  abund <- 100 * exp(-0.5 * ((mz - 200.0) / 0.002)^2)
  out   <- .Call("iso_centroid_Call", as.double(mz), as.double(abund),
                 as.integer(1), PACKAGE = "enviPat")
  expect_length(out[[1]], 1L)
  expect_equal(out[[1]][1], 200.0, tolerance = 1e-3)
})

test_that("iso_centroid_Call valley detected between two synthetic peaks", {
  mz2    <- seq(199.0, 202.0, by = 0.0001)
  abund2 <- 80  * exp(-0.5 * ((mz2 - 199.5) / 0.005)^2) +
            100 * exp(-0.5 * ((mz2 - 201.5) / 0.005)^2)
  out_v  <- .Call("iso_centroid_Call", as.double(mz2), as.double(abund2),
                  as.integer(2), PACKAGE = "enviPat")
  expect_gte(length(out_v[[1]]), 1L)
  expect_true(all(out_v[[1]] > 199.0 & out_v[[1]] < 202.0))
})

test_that("iso_centroid_Call flat-zero profile returns empty result", {
  mz_flat    <- seq(100.0, 101.0, by = 0.001)
  abund_flat <- rep(0.0, length(mz_flat))
  out_flat   <- .Call("iso_centroid_Call", as.double(mz_flat), as.double(abund_flat),
                      as.integer(0), PACKAGE = "enviPat")
  expect_length(out_flat[[1]], 0L)
})

# ── iso_profile_with_trace_Call ───────────────────────────────────────────────

test_that("iso_profile_with_trace_Call Gaussian returns list with m/z and abundance", {
  out_g <- .Call("iso_profile_with_trace_Call",
    as.integer(0), as.double(.peak_masses), as.double(.peak_abund),
    as.double(.trace), as.integer(1e5), as.double(0), as.integer(1),
    PACKAGE = "enviPat"
  )
  expect_type(out_g, "list")
  expect_length(out_g, 2L)
  expect_gt(length(out_g[[1]]), 0L)
  expect_true(all(out_g[[1]] >= 193.5 & out_g[[1]] <= 194.6))
  expect_lte(max(out_g[[2]]), 100.0 + 1e-6)
  expect_gt(max(out_g[[2]]), 0)
})

test_that("iso_profile_with_trace_Call Cauchy-Lorentz differs from Gaussian", {
  out_g <- .Call("iso_profile_with_trace_Call",
    as.integer(0), as.double(.peak_masses), as.double(.peak_abund),
    as.double(.trace), as.integer(1e5), as.double(0), as.integer(1),
    PACKAGE = "enviPat"
  )
  out_c <- .Call("iso_profile_with_trace_Call",
    as.integer(1), as.double(.peak_masses), as.double(.peak_abund),
    as.double(.trace), as.integer(1e5), as.double(0), as.integer(1),
    PACKAGE = "enviPat"
  )
  expect_gt(length(out_c[[1]]), 0L)
  expect_false(identical(out_g[[2]], out_c[[2]]))
})

test_that("iso_profile_with_trace_Call higher resolution gives narrower peak", {
  out_lo <- .Call("iso_profile_with_trace_Call",
    as.integer(0), as.double(.peak_masses), as.double(.peak_abund),
    as.double(.trace), as.integer(1e4), as.double(0), as.integer(1),
    PACKAGE = "enviPat"
  )
  out_hi <- .Call("iso_profile_with_trace_Call",
    as.integer(0), as.double(.peak_masses), as.double(.peak_abund),
    as.double(.trace), as.integer(1e6), as.double(0), as.integer(1),
    PACKAGE = "enviPat"
  )
  expect_gt(sum(out_lo[[2]] > 0.01), sum(out_hi[[2]] > 0.01))
})

# ── iso_pattern ───────────────────────────────────────────────────────────────

test_that("iso_pattern returns matrix with positive m/z and max abundance 100", {
  out <- .Call("iso_pattern",
    s1            = "C2H6O1",
    pl            = as.integer(1e6),
    t1            = as.double(0.001),
    iso_list_elem = iso_elem,
    iso_list_iso  = iso_iso,
    iso_list_mass = iso_mass,
    iso_list_abu  = iso_abu,
    rtm           = as.integer(0),
    rica          = as.integer(0),
    PACKAGE       = "enviPat"
  )
  expect_true(is.matrix(out))
  expect_gt(nrow(out), 0L)
  expect_true(all(out[, 1] > 0))
  expect_equal(max(out[, 2]), 100)
})

test_that("iso_pattern_4 returns same peak count and m/z values as iso_pattern", {
  args <- list(
    s1            = "C2H6O1",
    pl            = as.integer(1e6),
    t1            = as.double(0.001),
    iso_list_elem = iso_elem,
    iso_list_iso  = iso_iso,
    iso_list_mass = iso_mass,
    iso_list_abu  = iso_abu,
    rtm           = as.integer(0),
    rica          = as.integer(0),
    PACKAGE       = "enviPat"
  )
  out  <- do.call(.Call, c(list("iso_pattern"),  args))
  out4 <- do.call(.Call, c(list("iso_pattern_4"), args))
  expect_equal(nrow(out), nrow(out4))
  expect_equal(sort(out[, 1]), sort(out4[, 1]), tolerance = 1e-8)
})

test_that("iso_pattern returns empty first column for unknown formula", {
  out_bad <- .Call("iso_pattern",
    s1            = "ZZZZ999",
    pl            = as.integer(1e6),
    t1            = as.double(0.001),
    iso_list_elem = iso_elem,
    iso_list_iso  = iso_iso,
    iso_list_mass = iso_mass,
    iso_list_abu  = iso_abu,
    rtm           = as.integer(0),
    rica          = as.integer(0),
    PACKAGE       = "enviPat"
  )
  expect_equal(length(out_bad[[1]]), 0L)
})
