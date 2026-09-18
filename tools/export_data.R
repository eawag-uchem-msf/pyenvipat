#!/usr/bin/env Rscript
#
# Exports the four bundled data/*.rda datasets to plain CSV/TXT files that
# will be bundled into the Python port's package data (src/envipat/data/),
# see issue #6 of the R -> Python port (epic #2).
#
# Run from the repo root inside the dev container:
#   Rscript tools/export_data.R
#
# ---------------------------------------------------------------------------
# Precision-safe I/O
#
# write.csv() formats numerics via format(), which defaults to
# options("digits") (7 significant digits) and would silently truncate
# isotope masses that feed directly into the C core in later porting tasks.
# Every numeric value is instead written with 17 significant digits
# (the proven round-trip-exact bound for IEEE-754 doubles) via
# sprintf("%.17g", .), never via write.csv's own number formatting.
#
# This mirrors the approach already used by tools/make_golden.R (see issue
# #3 golden-fixture generator) for the same reason.
# ---------------------------------------------------------------------------

OUT <- "src/envipat/data"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

fmt_num <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.17g", x))

write_precise <- function(df, path) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    for (col in names(df)) {
        if (is.numeric(df[[col]])) df[[col]] <- fmt_num(df[[col]])
    }
    write.csv(df, path, row.names = FALSE)
}

cat("== enviPat data export (R -> CSV/TXT) ==\n")

## isotopes -> isotopes.csv ---------------------------------------------------
## columns: element, isotope, mass, abundance, ratioC (unchanged from the
## .rda; mass/abundance are doubles and get full-precision export).
load("data/isotopes.rda")
out_isotopes <- file.path(OUT, "isotopes.csv")
write_precise(isotopes, out_isotopes)
cat(sprintf("wrote %s (%d rows)\n", out_isotopes, nrow(isotopes)))

## chemforms -> chemforms.txt -------------------------------------------------
## one formula per line, no header, no quoting needed (formulas contain no
## commas or newlines).
load("data/chemforms.rda")
out_chemforms <- file.path(OUT, "chemforms.txt")
writeLines(chemforms, out_chemforms)
cat(sprintf("wrote %s (%d lines)\n", out_chemforms, length(chemforms)))

## adducts -> adducts.csv -----------------------------------------------------
## 9 columns, same names/order as R/data-adducts.R:
## Name, calc, Charge, Mult, Mass, Ion_mode, Formula_add, Formula_ded, Multi
load("data/adducts.rda")
out_adducts <- file.path(OUT, "adducts.csv")
write_precise(adducts, out_adducts)
cat(sprintf("wrote %s (%d rows)\n", out_adducts, nrow(adducts)))

## resolution_list -> resolution_list.csv ------------------------------------
## Long format: one row per curve point, `name` identifies which curve
## (instrument/resolution setting) the (mz, R) pair belongs to. The original
## per-curve data frames use the literal column name "m/z", which is legal
## but awkward as a CSV/pandas column name, so it is renamed to "mz" on
## export (documented here and in src/envipat/data.py).
##
## Note: not all elements of resolution_list are data.frames - a handful are
## plain numeric matrices with the same two columns ("m/z", "R"). `df[[1]]`
## on a matrix returns a single flattened element (not a column!), so
## columns are extracted with `df[, 1]` / `df[, 2]`, which correctly returns
## the column vector for both matrices and data.frames.
load("data/resolution_list.rda")
curve_names <- names(resolution_list)
long <- do.call(rbind, lapply(curve_names, function(nm) {
    df <- resolution_list[[nm]]
    data.frame(
        name = nm,
        mz = df[, 1],   # originally "m/z"
        R = df[, 2],
        stringsAsFactors = FALSE
    )
}))
out_resolution <- file.path(OUT, "resolution_list.csv")
write_precise(long, out_resolution)
cat(sprintf(
    "wrote %s (%d rows, %d curves)\n",
    out_resolution, nrow(long), length(curve_names)
))

cat("done.\n")
