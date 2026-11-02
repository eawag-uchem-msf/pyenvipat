#!/usr/bin/env Rscript
#
# Round-trip verification for tools/export_data.R (issue #6): reloads the
# CSV/TXT files written into src/envipat/data/ and checks them against the
# original data/*.rda objects with identical() (or an explicitly documented
# equivalent exact check, where identical()'s notion of "identical" doesn't
# survive a text format - see the resolution_list note below).
#
# Run from the repo root inside the dev container, after tools/export_data.R:
#   Rscript tools/verify_export.R

DATA <- "src/envipat/data"
ok_all <- TRUE
check <- function(label, ok) {
    cat(sprintf("[%s] %s\n", if (ok) "OK" else "FAIL", label))
    ok_all <<- ok_all && ok
}

## isotopes -------------------------------------------------------------------
## Note: colClasses is deliberately NOT used here - passing colClasses =
## "numeric"/"integer" makes read.table() hand quoted fields (our
## write.csv(quote = TRUE) output) straight to a low-level scan() that does
## not strip quotes, and errors out. Auto-detection (the default) tokenizes
## quotes correctly and then type-converts, so columns are read generically
## and then cast to match the original's classes explicitly.
load("data/isotopes.rda")
reload_isotopes <- read.csv(file.path(DATA, "isotopes.csv"), stringsAsFactors = FALSE)
reload_isotopes$ratioC <- as.integer(reload_isotopes$ratioC)
check("isotopes identical() to original .rda", identical(isotopes, reload_isotopes))

## chemforms --------------------------------------------------------------------
load("data/chemforms.rda")
reload_chemforms <- readLines(file.path(DATA, "chemforms.txt"))
check("chemforms identical() to original .rda", identical(chemforms, reload_chemforms))

## adducts ----------------------------------------------------------------------
load("data/adducts.rda")
reload_adducts <- read.csv(file.path(DATA, "adducts.csv"), stringsAsFactors = FALSE)
# Charge/Mult/Multi happen to hold only whole numbers, so read.csv's
# auto-detection (type.convert) infers class integer for them, even though
# the original columns are class numeric (double) - again a distinction
# with no CSV representation (see the resolution_list note above). Cast
# back to numeric to compare on value, not on this incidental R-side type.
for (col in c("Charge", "Mult", "Multi")) {
    reload_adducts[[col]] <- as.numeric(reload_adducts[[col]])
}
# The original adducts.rda also carries non-sequential character row.names
# (e.g. "46","47",...,"11","7",... - a leftover from how the object was
# assembled upstream, unrelated to its 9 documented columns). Row NAMES
# are not part of the CSV/9-column deliverable and CSV has no row-name
# concept at all, so row.names are reset to the default on both sides
# before comparing; row ORDER (which does matter) is preserved and still
# checked, since it drives every other column's identical() comparison.
adducts_cmp <- adducts
row.names(adducts_cmp) <- NULL
row.names(reload_adducts) <- NULL
check("adducts identical() to original .rda (row order + all 9 columns; row.names labels intentionally not preserved)",
      identical(adducts_cmp, reload_adducts))

## resolution_list ----------------------------------------------------------------
## Not all 31 original list elements are the same R container: most are
## data.frames, but a handful (Sciex_*/Agilent_* curves) are plain numeric
## matrices - a distinction with no representation in a CSV file, just like
## R's integer-vs-double distinction for whole-number columns (e.g. many
## curves' "R" column is class integer, others class numeric, purely
## depending on how each curve was originally entered - not a meaningful
## semantic difference for m/z or resolution values). identical() on the
## containers themselves would therefore spuriously fail even on a perfect
## export, so this check instead compares column VALUES exactly (via
## as.numeric() coercion on both sides, using identical() for the
## comparison itself, i.e. full bit-exact equality, no tolerance) per curve,
## in original curve order - the same "values, not container/dtype" standard
## the issue explicitly asks for on the Python side.
load("data/resolution_list.rda")
long <- read.csv(file.path(DATA, "resolution_list.csv"), stringsAsFactors = FALSE)

resolution_ok <- TRUE
curve_names <- names(resolution_list)
if (!identical(sort(unique(long$name)), sort(curve_names))) {
    resolution_ok <- FALSE
    cat("  curve name set mismatch\n")
}
for (nm in curve_names) {
    orig <- resolution_list[[nm]]
    sub <- long[long$name == nm, c("mz", "R")]
    orig_mz <- as.numeric(orig[, 1])
    orig_r <- as.numeric(orig[, 2])
    ok_nrow <- nrow(sub) == nrow(orig)
    ok_mz <- ok_nrow && identical(orig_mz, as.numeric(sub$mz))
    ok_r <- ok_nrow && identical(orig_r, as.numeric(sub$R))
    if (!(ok_nrow && ok_mz && ok_r)) {
        resolution_ok <- FALSE
        cat(sprintf(
            "  mismatch in curve %s: nrow_ok=%s mz_ok=%s R_ok=%s\n",
            nm, ok_nrow, ok_mz, ok_r
        ))
    }
}
check("resolution_list: all 31 curves match original values exactly", resolution_ok)

cat("\n")
if (ok_all) {
    cat("All round-trip checks passed.\n")
} else {
    cat("Some round-trip checks FAILED.\n")
    quit(status = 1)
}
