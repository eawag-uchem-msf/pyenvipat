#!/usr/bin/env Rscript
#
# Dumps the installed enviPat package's `isotopes` data(...) table to a CSV
# at full double precision (%.17g, round-trip exact for IEEE-754 doubles --
# same convention tools/make_golden.R uses), for the Python C extension's
# own golden-fixture verification (see tests/test_core.py).
#
# Several tests/golden/c_direct/ cases (issue #3) call iso_pattern/
# iso_pattern_4 with the *default* `isotopes` table rather than a
# committed, reduced one (only tests/golden/c_direct/case_003_input_
# isotopes.csv, H/C/O only, is committed as a fixture input). Since issue
# #5 branches off #4, not #6 (which exports data/isotopes.rda to CSV), this
# script is the only source of the full table available to it: it re-reads
# the same installed R package/data(isotopes) the golden fixtures themselves
# were generated from, rather than re-deriving or hand-copying isotope mass/
# abundance values (transcription of ~300 rows would itself be a source of
# error a "bit-exact" claim can't tolerate).
#
# Usage: Rscript tools/dump_isotopes_precise.R <output.csv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("usage: Rscript tools/dump_isotopes_precise.R <output.csv>")
}

suppressMessages(library(enviPat))
data(isotopes)

fmt_num <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.17g", x))

df <- as.data.frame(isotopes, stringsAsFactors = FALSE)
for (col in names(df)) {
    if (is.numeric(df[[col]])) df[[col]] <- fmt_num(df[[col]])
}
write.csv(df, args[1], row.names = FALSE)
