#!/usr/bin/env Rscript
#
# Full end-to-end pipeline dump for issue #10's parity verification.
#
# Unlike tools/dump_isowrap_resolution.R (issue #8), which only needed a
# resolution *vector* to feed into Python's isowrap(resolution=...)
# (issue #9's getR port wasn't merged into that branch), this script
# exercises R's real isowrap() end-to-end -- letting it call getR()
# *internally* from a resmass calibration table (nknots=13, spar=0.1,
# R's own getR() defaults) -- across the full cross-product issue #10
# asks for: the full data(chemforms) set (107 formulas, confirmed
# byte-identical to Python's envipat.data.chemforms -- see MR !9's
# post-review update) x all 31 resolution_list curves x algo {1, 2} x
# detect {centroid, intensoid, valley}. An earlier revision of this
# script used a hand-picked 8-formula subset; the reviewer asked for
# full coverage instead (MR !9 post-review), and all 107 formulas pass
# check_chemform() cleanly (verified below), so no formulas end up
# skipped in practice.
#
# The Python side (tests/python/test_end_to_end.py) does NOT replay
# these exact isowrap() calls -- it instead calls envipat.resolution.
# get_resolution() itself (issue #9's from-scratch getR port) to derive
# its own resolution per curve, then composes isopattern -> envelope ->
# vdetect independently. That is the point of this dump: both sides
# derive their own resolution via their own getR/get_resolution
# implementation, so comparing final output is a genuine end-to-end
# check of the whole chain, not just of isowrap's own composition (which
# issue #8 already verified with a pre-computed resolution vector).
#
# Output: <output_dir>/resolution.csv (curve, formula, resolution) and
# <output_dir>/results.csv (curve, algo, detect, formula, row, mz,
# abundance, is_error) -- long format, one row per output peak, matching
# the m/z-then-abundance column order every step in this pipeline uses.
#
# Usage: Rscript tools/dump_full_pipeline.R <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("usage: Rscript tools/dump_full_pipeline.R <output_dir>")
}
outdir <- args[1]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(enviPat))
data(isotopes)
data(chemforms)
data(resolution_list)

fmt_num <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.17g", x))

# Full data(chemforms) set (107 formulas as of this package version).
# check_chemform() failures are skipped explicitly (logged to stderr and
# dropped from the matrix) rather than aborting the whole run over one
# bad formula -- as of this revision all 107 pass, but this keeps the
# script robust if that data set ever changes.
formulas <- as.character(chemforms)

checked <- check_chemform(isotopes, formulas)
bad <- checked[, 1]
if (any(bad)) {
    message(
        "skipping ", sum(bad), " formula(s) that failed check_chemform: ",
        paste(formulas[bad], collapse = ", ")
    )
}
formulas <- formulas[!bad]
new_formulas <- checked[!bad, 2]
# getR() itself stops on any(checked[, 1] == TRUE), so it needs the
# filtered table too, not just the filtered formula names.
checked <- checked[!bad, ]

# isopattern doesn't depend on the resolution curve -- compute both algos'
# patterns once, outside the curve loop, and reuse them for every curve.
patterns <- list()
for (algo in c(1, 2)) {
    patterns[[as.character(algo)]] <- isopattern(
        isotopes, new_formulas,
        threshold = 0.1, charge = 1, algo = algo, verbose = FALSE
    )
}

res_rows <- list()
result_rows <- list()
detect_modes <- c("centroid", "intensoid", "valley")

for (curve in names(resolution_list)) {
    resmass <- resolution_list[[curve]]
    # getR() (both R's and get_resolution's Python port) rejects
    # nknots > nrow(resmass) -- 2 of the 31 curves
    # (Agilent_QTOF6550_low_extended2GHz_highRes,
    # Agilent_QTOF6550_low_highRes4GHz_highRes) only have 12 calibration
    # rows, below getR's own nknots=13 default. This is a shared
    # constraint of the algorithm itself (identical check in both
    # languages), not an R/Python divergence, so both sides here use
    # nknots=min(13, nrow(resmass)) rather than skipping those curves.
    nknots <- min(13, nrow(resmass))
    resolution <- getR(checked, resmass = resmass, nknots = nknots, spar = 0.1, plotit = FALSE)

    for (i in seq_along(new_formulas)) {
        res_rows[[length(res_rows) + 1]] <- data.frame(
            curve = curve, formula = new_formulas[i], nknots = nknots,
            resolution = fmt_num(resolution[i]),
            stringsAsFactors = FALSE
        )
    }

    for (algo in c(1, 2)) {
        pattern <- patterns[[as.character(algo)]]
        profiles <- envelope(
            pattern,
            ppm = FALSE, dmz = "get", frac = 0.25, env = "Gaussian",
            resolution = resolution, verbose = FALSE
        )
        for (detect in detect_modes) {
            centro <- vdetect(profiles, detect = detect, verbose = FALSE)
            for (i in seq_along(new_formulas)) {
                fm <- new_formulas[i]
                df <- centro[[fm]]
                if (is.null(df) || identical(df, "error") || nrow(df) == 0) {
                    result_rows[[length(result_rows) + 1]] <- data.frame(
                        curve = curve, algo = algo, detect = detect, formula = fm,
                        row = NA_integer_, mz = NA_character_, abundance = NA_character_,
                        is_error = TRUE, stringsAsFactors = FALSE
                    )
                } else {
                    for (r in seq_len(nrow(df))) {
                        result_rows[[length(result_rows) + 1]] <- data.frame(
                            curve = curve, algo = algo, detect = detect, formula = fm,
                            row = r, mz = fmt_num(df[r, 1]), abundance = fmt_num(df[r, 2]),
                            is_error = FALSE, stringsAsFactors = FALSE
                        )
                    }
                }
            }
        }
    }
}

write.csv(do.call(rbind, res_rows), file.path(outdir, "resolution.csv"), row.names = FALSE)
write.csv(do.call(rbind, result_rows), file.path(outdir, "results.csv"), row.names = FALSE)
cat("wrote", nrow(do.call(rbind, res_rows)), "resolution rows and",
    nrow(do.call(rbind, result_rows)), "result rows to", outdir, "\n")
