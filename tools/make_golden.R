#!/usr/bin/env Rscript
#
# Generates golden reference fixtures from the current, unmodified R package
# (issue #3 of the R -> Python port, see design/python-port-plan.md). These
# fixtures define "correct" for every later porting task and, once R is
# removed, become the permanent Python regression suite.
#
# Run from the repo root inside the dev container:
#   R CMD INSTALL .
#   Rscript tools/make_golden.R
#
# Determinism: nothing here is randomized. Re-running must reproduce
# byte-identical output (verified by CI / manual `git status` check).

suppressMessages(library(enviPat))
data(isotopes); data(chemforms); data(adducts); data(resolution_list)

OUT <- "tests/golden"
unlink(OUT, recursive = TRUE)
dir.create(OUT, recursive = TRUE)

## ---------------------------------------------------------------------------
## Precision-safe I/O helpers
##
## write.csv() formats numerics via options("digits") (default 7), which
## would silently truncate isotope masses that feed straight into the C core.
## Every numeric is written with 17 significant digits (round-trip exact for
## IEEE-754 doubles) via sprintf, never via write.csv's own formatting.
## ---------------------------------------------------------------------------

fmt_num <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.17g", x))

write_precise <- function(df, path) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    for (col in names(df)) {
        if (is.numeric(df[[col]])) df[[col]] <- fmt_num(df[[col]])
    }
    write.csv(df, path, row.names = FALSE)
}

write_vec <- function(x, path, colname = "value") {
    d <- data.frame(x, stringsAsFactors = FALSE)
    names(d) <- colname
    write_precise(d, path)
}

rbind_fill <- function(df_list) {
    # Base-R rbind that unions columns and fills missing ones with NA -
    # isopattern's per-formula matrices have different isotope columns
    # depending on which elements/isotopes are actually present.
    all_cols <- unique(unlist(lapply(df_list, names)))
    df_list2 <- lapply(df_list, function(df) {
        missing <- setdiff(all_cols, names(df))
        for (m in missing) df[[m]] <- NA
        df[all_cols]
    })
    do.call(rbind, df_list2)
}

params_str <- function(...) {
    args <- list(...)
    paste(sprintf("%s=%s", names(args), vapply(args, as.character, "")), collapse = ";")
}

manifests <- new.env()
add_case <- function(section, case_id, fn, params, output_file, input_file = "", notes = "") {
    row <- data.frame(
        case_id = case_id, fn = fn, params = params,
        output_file = output_file, input_file = input_file, notes = notes,
        stringsAsFactors = FALSE
    )
    if (is.null(manifests[[section]])) manifests[[section]] <- list()
    manifests[[section]][[length(manifests[[section]]) + 1]] <- row
}
flush_manifest <- function(section, dir) {
    df <- do.call(rbind, manifests[[section]])
    write.csv(df, file.path(dir, "cases.csv"), row.names = FALSE)
}

stack_pattern_list <- function(lst) {
    parts <- lapply(names(lst), function(nm) {
        el <- lst[[nm]]
        if (identical(el, "error")) {
            data.frame(compound = nm, error = TRUE, stringsAsFactors = FALSE)
        } else {
            df <- as.data.frame(el, stringsAsFactors = FALSE)
            cbind(compound = nm, error = FALSE, df, stringsAsFactors = FALSE)
        }
    })
    rbind_fill(parts)
}

cat("== enviPat golden fixture generator ==\n")

## =============================================================================
## 1. Direct .Call entry points (bypass R wrappers entirely)
## =============================================================================
dir_c <- file.path(OUT, "c_direct"); dir.create(dir_c, recursive = TRUE)
n <- 0

call_iso_pattern <- function(entry, s1, t1, rtm, rica, iso = isotopes) {
    .Call(entry, s1 = as.character(s1), pl = as.integer(1e6), t1 = as.double(t1),
        iso_list_elem = as.character(iso[, 1]), iso_list_iso = as.character(iso[, 2]),
        iso_list_mass = as.numeric(iso[, 3]), iso_list_abu = as.numeric(iso[, 4]),
        rtm = as.integer(rtm), rica = as.integer(rica), PACKAGE = "enviPat")
}

# 1a. iso_pattern / iso_pattern_4, normal formula, both algorithms
for (entry in c("iso_pattern", "iso_pattern_4")) {
    out <- call_iso_pattern(entry, "C2H4", 0.1, 0, 0)
    n <- n + 1; id <- sprintf("case_%03d", n)
    write_precise(out, file.path(dir_c, paste0(id, ".csv")))
    add_case("c_direct", id, entry, params_str(s1 = "C2H4", pl = 1e6, t1 = 0.1, rtm = 0, rica = 0),
        paste0(id, ".csv"))
}

# 1b. custom/reduced isotope table (only H, C, O rows) - proves the C code
# honours the table argument rather than a hardcoded one.
iso_reduced <- isotopes[isotopes[, 1] %in% c("H", "C", "O"), ]
out <- call_iso_pattern("iso_pattern", "C2H4O1", 0.1, 0, 0, iso = iso_reduced)
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(out, file.path(dir_c, paste0(id, ".csv")))
write_precise(iso_reduced, file.path(dir_c, paste0(id, "_input_isotopes.csv")))
add_case("c_direct", id, "iso_pattern", params_str(s1 = "C2H4O1", pl = 1e6, t1 = 0.1, rtm = 0, rica = 0),
    paste0(id, ".csv"), input_file = paste0(id, "_input_isotopes.csv"), notes = "reduced isotope table (H,C,O only)")

# 1c. rica=1 (return_iso_calc_amount) path - returns a named int vector, not a matrix
out <- call_iso_pattern("iso_pattern", "C2H4", 0.1, 0, 1)
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(data.frame(element = names(out), count = as.integer(out)), file.path(dir_c, paste0(id, ".csv")))
add_case("c_direct", id, "iso_pattern", params_str(s1 = "C2H4", pl = 1e6, t1 = 0.1, rtm = 0, rica = 1),
    paste0(id, ".csv"), notes = "rica=1: returns named int vector of element counts, not isotope matrix")

# 1d. invalid formula -> NULL (out[[1]] has length 0 -> R wrapper reports "error")
out <- call_iso_pattern("iso_pattern", "ZZZ99", 0.1, 0, 0)
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(data.frame(is_null = is.null(out)), file.path(dir_c, paste0(id, ".csv")))
add_case("c_direct", id, "iso_pattern", params_str(s1 = "ZZZ99", pl = 1e6, t1 = 0.1, rtm = 0, rica = 0),
    paste0(id, ".csv"), notes = "invalid formula -> C returns NULL")

# 1e. rtm (rel_to) sweep on a formula with several isotope peaks
for (rtm in 0:4) {
    out <- call_iso_pattern("iso_pattern", "C10H15Cl2", 0.001, rtm, 0)
    n <- n + 1; id <- sprintf("case_%03d", n)
    write_precise(out, file.path(dir_c, paste0(id, ".csv")))
    add_case("c_direct", id, "iso_pattern", params_str(s1 = "C10H15Cl2", pl = 1e6, t1 = 0.001, rtm = rtm, rica = 0),
        paste0(id, ".csv"))
}

# 1f. iso_ppm_Call - a few start/end/ppm combinations
ppm_cases <- list(
    c(100.0, 101.0, 5.0), c(100.0, 101.0, 1.0), c(100.0, 110.0, 5.0), c(100.0, 105.0, 3.0)
)
for (pc in ppm_cases) {
    out <- .Call("iso_ppm_Call", as.double(pc[1]), as.double(pc[2]), as.double(pc[3]), PACKAGE = "enviPat")
    n <- n + 1; id <- sprintf("case_%03d", n)
    write_vec(out, file.path(dir_c, paste0(id, ".csv")), colname = "mz")
    add_case("c_direct", id, "iso_ppm_Call", params_str(start = pc[1], end = pc[2], ppm = pc[3]),
        paste0(id, ".csv"))
}

# 1g. iso_profile_with_trace_Call - Gaussian/CauchyLorentz, filter on/off, t1 sweep
trace1 <- seq(193.5, 194.6, by = 0.001)
profile_grid <- expand.grid(type1 = c(0, 1), filter = c(0, 1), t1 = c(0, 0.01), stringsAsFactors = FALSE)
profile_outputs <- list()
for (i in seq_len(nrow(profile_grid))) {
    g <- profile_grid[i, ]
    out <- .Call("iso_profile_with_trace_Call",
        type1 = as.integer(g$type1), f1 = as.double(194.08), a1 = as.double(100.0),
        tr1 = as.double(trace1), r1 = as.integer(1e5), t1 = as.double(g$t1),
        filter = as.integer(g$filter), PACKAGE = "enviPat")
    n <- n + 1; id <- sprintf("case_%03d", n)
    write_precise(data.frame(mz = out[[1]], abundance = out[[2]]), file.path(dir_c, paste0(id, ".csv")))
    write_vec(trace1, file.path(dir_c, paste0(id, "_input_trace.csv")), colname = "mz")
    add_case("c_direct", id, "iso_profile_with_trace_Call",
        params_str(type1 = g$type1, f1 = 194.08, a1 = 100.0, r1 = 1e5, t1 = g$t1, filter = g$filter),
        paste0(id, ".csv"), input_file = paste0(id, "_input_trace.csv"))
    profile_outputs[[paste(g$type1, g$filter, g$t1)]] <- out
}

# 1h. iso_centroid_Call - all 3 detect types, on the filter=1/t1=0/Gaussian profile above
base_profile <- profile_outputs[["0 1 0"]]
for (t1 in 0:2) {
    out <- .Call("iso_centroid_Call", f1 = as.double(base_profile[[1]]), a1 = as.double(base_profile[[2]]),
        t1 = as.integer(t1), PACKAGE = "enviPat")
    n <- n + 1; id <- sprintf("case_%03d", n)
    write_precise(data.frame(mz = out[[1]], abundance = out[[2]]), file.path(dir_c, paste0(id, ".csv")))
    add_case("c_direct", id, "iso_centroid_Call", params_str(t1 = t1, profile = "gaussian_filter1_t0"),
        paste0(id, ".csv"), input_file = "case_017.csv")
}

# 1i. iso_centroid_Call on a profile with a real valley and tied plateau
# intensities, to pin tie-breaking behaviour and the documented dead-code
# read (c_a[j] read immediately after j++, before being written for that
# index - see design/python-port-plan.md quirks list).
tied_out <- .Call("iso_profile_with_trace_Call", type1 = as.integer(0),
    f1 = as.double(c(194.0, 194.2)), a1 = as.double(c(100.0, 100.0)),
    tr1 = as.double(seq(193.5, 194.7, by = 0.001)), r1 = as.integer(2e4),
    t1 = as.double(0), filter = as.integer(1), PACKAGE = "enviPat")
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(data.frame(mz = tied_out[[1]], abundance = tied_out[[2]]), file.path(dir_c, paste0(id, ".csv")))
tied_id <- id
for (t1 in 0:2) {
    out <- .Call("iso_centroid_Call", f1 = as.double(tied_out[[1]]), a1 = as.double(tied_out[[2]]),
        t1 = as.integer(t1), PACKAGE = "enviPat")
    n <- n + 1; id <- sprintf("case_%03d", n)
    write_precise(data.frame(mz = out[[1]], abundance = out[[2]]), file.path(dir_c, paste0(id, ".csv")))
    add_case("c_direct", id, "iso_centroid_Call", params_str(t1 = t1, profile = "two_equal_peaks_valley"),
        paste0(id, ".csv"), input_file = paste0(tied_id, ".csv"))
}

flush_manifest("c_direct", dir_c)
cat(sprintf("c_direct: %d cases\n", n))

## =============================================================================
## Shared formula pools
## =============================================================================
ck_all <- check_chemform(isotopes, chemforms)
mass_order <- order(ck_all[, 3])

pool10 <- unique(c(
    chemforms[mass_order[1]],                                   # min mass
    chemforms[mass_order[length(mass_order)]],                  # max mass
    chemforms[mass_order[round(length(mass_order) / 2)]],        # median mass
    chemforms[mass_order[round(length(mass_order) / 4)]],        # another mid-mass
    (chemforms[grepl("Cl[2-9]", chemforms)])[1],                 # Cl-heavy
    (chemforms[grepl("Br[1-9]", chemforms)])[1],                 # Br-heavy
    "C1H4",                                                       # minimal formula
    "C50H80N10O15S2",                                             # synthetic large/peptide-like
    "[13]C2C8D2H8[15]N2N2[18]O2O3",                               # isotope-labeled (mixed D/heavy/light)
    "C20Cl10"                                                     # many isotopes of one element
))
pool3 <- pool10[c(1, 7, 9)]  # small, minimal, isotope-labeled - keep envelope/vdetect output size sane

## =============================================================================
## 2. isopattern
## =============================================================================
dir_ip <- file.path(OUT, "isopattern"); dir.create(dir_ip, recursive = TRUE)
n <- 0

# 2a. Grid A: pool10 x algo{1,2} x rel_to{0..4} x threshold{0.1, 0.001}, charge=FALSE
grid_a <- expand.grid(algo = c(1, 2), rel_to = 0:4, threshold = c(0.1, 0.001))
for (i in seq_len(nrow(grid_a))) {
    g <- grid_a[i, ]
    pat <- isopattern(isotopes, pool10, threshold = g$threshold, charge = FALSE,
        algo = g$algo, rel_to = g$rel_to, verbose = FALSE)
    n <- n + 1; id <- sprintf("case_%03d", n)
    write_precise(stack_pattern_list(pat), file.path(dir_ip, paste0(id, ".csv")))
    add_case("isopattern", id, "isopattern",
        params_str(algo = g$algo, rel_to = g$rel_to, threshold = g$threshold, charge = "FALSE", n_formulas = length(pool10)),
        paste0(id, ".csv"), notes = "grid_a")
}

# 2b. Grid B: pool10 x charge {1, 2, -1, per-formula vector}, algo=1 rel_to=0 threshold=0.1
charge_cases <- list("1" = 1, "2" = 2, "-1" = -1)
for (nm in names(charge_cases)) {
    pat <- isopattern(isotopes, pool10, threshold = 0.1, charge = charge_cases[[nm]],
        algo = 1, rel_to = 0, verbose = FALSE)
    n <- n + 1; id <- sprintf("case_%03d", n)
    write_precise(stack_pattern_list(pat), file.path(dir_ip, paste0(id, ".csv")))
    add_case("isopattern", id, "isopattern",
        params_str(algo = 1, rel_to = 0, threshold = 0.1, charge = nm, n_formulas = length(pool10)),
        paste0(id, ".csv"), notes = "grid_b_scalar_charge")
}
charge_vec <- rep(c(1, 2, -1), length.out = length(pool10))
pat <- isopattern(isotopes, pool10, threshold = 0.1, charge = charge_vec, algo = 1, rel_to = 0, verbose = FALSE)
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(stack_pattern_list(pat), file.path(dir_ip, paste0(id, ".csv")))
add_case("isopattern", id, "isopattern",
    params_str(algo = 1, rel_to = 0, threshold = 0.1, charge = paste(charge_vec, collapse = "/"), n_formulas = length(pool10)),
    paste0(id, ".csv"), notes = "grid_b_per_formula_charge_vector")

# 2c. Broad pass: all 107 chemforms, default algo/rel_to, two thresholds
for (thr in c(0.1, 0.001)) {
    pat <- isopattern(isotopes, chemforms, threshold = thr, charge = FALSE, algo = 1, rel_to = 0, verbose = FALSE)
    n <- n + 1; id <- sprintf("case_%03d", n)
    write_precise(stack_pattern_list(pat), file.path(dir_ip, paste0(id, ".csv")))
    add_case("isopattern", id, "isopattern",
        params_str(algo = 1, rel_to = 0, threshold = thr, charge = "FALSE", n_formulas = length(chemforms)),
        paste0(id, ".csv"), notes = "broad_pass_all_chemforms")
}

# 2d. threshold=0: full, unpruned fine structure on two small formulas
pat <- isopattern(isotopes, c("C1H4", pool10[1]), threshold = 0, charge = FALSE, algo = 1, rel_to = 0, verbose = FALSE)
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(stack_pattern_list(pat), file.path(dir_ip, paste0(id, ".csv")))
add_case("isopattern", id, "isopattern",
    params_str(algo = 1, rel_to = 0, threshold = 0, charge = "FALSE", n_formulas = 2),
    paste0(id, ".csv"), notes = "threshold_0_full_fine_structure")

flush_manifest("isopattern", dir_ip)
cat(sprintf("isopattern: %d cases\n", n))

## =============================================================================
## 3. envelope + vdetect
## =============================================================================
dir_env <- file.path(OUT, "envelope"); dir.create(dir_env, recursive = TRUE)
dir_vd  <- file.path(OUT, "vdetect");  dir.create(dir_vd,  recursive = TRUE)
n_env <- 0; n_vd <- 0

pattern3 <- isopattern(isotopes, pool3, threshold = 0.1, charge = FALSE, algo = 1, verbose = FALSE)

dmz_modes <- list(
    list(name = "get", ppm = FALSE, dmz = "get"),
    list(name = "abs", ppm = FALSE, dmz = 0.001),
    list(name = "ppm", ppm = TRUE,  dmz = 5)
)
env_grid <- expand.grid(env = c("Gaussian", "CauchyLorentz"), resolution = c(2e4, 4e4), dmz_i = 1:3, stringsAsFactors = FALSE)

for (i in seq_len(nrow(env_grid))) {
    g <- env_grid[i, ]
    dm <- dmz_modes[[g$dmz_i]]
    prof <- envelope(pattern3, ppm = dm$ppm, dmz = dm$dmz, frac = 1/4, env = g$env,
        resolution = g$resolution, verbose = FALSE)
    n_env <- n_env + 1; eid <- sprintf("case_%03d", n_env)
    write_precise(stack_pattern_list(prof), file.path(dir_env, paste0(eid, ".csv")))
    add_case("envelope", eid, "envelope",
        params_str(env = g$env, resolution = g$resolution, ppm = dm$ppm, dmz = dm$name, frac = "0.25", n_formulas = length(pool3)),
        paste0(eid, ".csv"))

    for (detect in c("centroid", "intensoid", "valley")) {
        det <- vdetect(prof, detect = detect, plotit = FALSE, verbose = FALSE)
        n_vd <- n_vd + 1; vid <- sprintf("case_%03d", n_vd)
        write_precise(stack_pattern_list(det), file.path(dir_vd, paste0(vid, ".csv")))
        add_case("vdetect", vid, "vdetect", params_str(detect = detect, source_envelope_case = eid),
            paste0(vid, ".csv"), input_file = file.path("..", "envelope", paste0(eid, ".csv")))
    }
}

flush_manifest("envelope", dir_env)
flush_manifest("vdetect", dir_vd)
cat(sprintf("envelope: %d cases, vdetect: %d cases\n", n_env, n_vd))

## =============================================================================
## 4. getR - one file per resolution curve, one column per (nknots, spar) combo
## =============================================================================
dir_gr <- file.path(OUT, "getR"); dir.create(dir_gr, recursive = TRUE)
n <- 0
combos <- list(c(13, 0.1), c(6, 0.2), c(13, 0.5), c(20, 0.1))

for (curve_name in names(resolution_list)) {
    resmass <- resolution_list[[curve_name]]
    mass_grid <- seq(min(resmass[, 1]), max(resmass[, 1]), length.out = 15)
    checked <- data.frame(warning = FALSE, formula = NA_character_, mass = mass_grid)

    result <- data.frame(mass = mass_grid)
    used_combos <- c()
    for (combo in combos) {
        nknots <- combo[1]; spar <- combo[2]
        if (nknots > nrow(resmass)) next
        res <- getR(checked, resmass, nknots = nknots, spar = spar, plotit = FALSE)
        colname <- sprintf("resolution_nknots%d_spar%s", nknots, spar)
        result[[colname]] <- res
        used_combos <- c(used_combos, sprintf("(%d,%s)", nknots, spar))
    }
    n <- n + 1; id <- sprintf("case_%03d", n)
    write_precise(result, file.path(dir_gr, paste0(id, ".csv")))
    add_case("getR", id, "getR",
        params_str(curve = curve_name, n_points = nrow(resmass), combos = paste(used_combos, collapse = "/")),
        paste0(id, ".csv"))
}
flush_manifest("getR", dir_gr)
cat(sprintf("getR: %d cases\n", n))

## =============================================================================
## 5. isowrap - end-to-end pipeline reference
## =============================================================================
dir_iw <- file.path(OUT, "isowrap"); dir.create(dir_iw, recursive = TRUE)
n <- 0

pool5 <- pool10[1:5]
checked5 <- check_chemform(isotopes, pool5)
curve_names_3 <- c("Elite_R240000@400", "Q-Exactive,ExactivePlus_R70000@200", "QTOF_XevoG2-S_R25000@200")

iw_grid <- expand.grid(curve = curve_names_3, algo = c(1, 2), detect = c("centroid", "intensoid", "valley"), stringsAsFactors = FALSE)
for (i in seq_len(nrow(iw_grid))) {
    g <- iw_grid[i, ]
    resmass <- resolution_list[[g$curve]]
    out <- isowrap(isotopes, checked5, resmass, resolution = FALSE, nknots = 6, spar = 0.2,
        threshold = 0.1, charge = 1, algo = g$algo, ppm = FALSE, dmz = "get", frac = 1/4,
        env = "Gaussian", detect = g$detect, plotit = FALSE, verbose = FALSE)
    n <- n + 1; id <- sprintf("case_%03d", n)
    write_precise(stack_pattern_list(out), file.path(dir_iw, paste0(id, ".csv")))
    add_case("isowrap", id, "isowrap",
        params_str(curve = g$curve, algo = g$algo, detect = g$detect, nknots = 6, spar = 0.2, threshold = 0.1, charge = 1),
        paste0(id, ".csv"))
}
flush_manifest("isowrap", dir_iw)
cat(sprintf("isowrap: %d cases\n", n))

## =============================================================================
## 6. check_several
## =============================================================================
dir_cs <- file.path(OUT, "check_several"); dir.create(dir_cs, recursive = TRUE)
n <- 0

pat_overlap    <- isopattern(isotopes, c("C10H10", "C9H10[13]C1"), threshold = 0.1, verbose = FALSE)
pat_no_overlap <- isopattern(isotopes, pool3, threshold = 0.1, verbose = FALSE)

cs_cases <- list(
    list(name = "ppmFALSE_overlap",    pattern = pat_overlap,    dmz = 0.01, ppm = FALSE),
    list(name = "ppmFALSE_no_overlap", pattern = pat_no_overlap, dmz = 1e-6, ppm = FALSE),
    list(name = "ppmTRUE_no_overlap",  pattern = pat_no_overlap, dmz = 1e-3, ppm = TRUE)
)
for (c1 in cs_cases) {
    res <- check_several(c1$pattern, dmz = c1$dmz, ppm = c1$ppm)
    n <- n + 1; id <- sprintf("case_%03d", n)
    write_precise(res, file.path(dir_cs, paste0(id, ".csv")))
    add_case("check_several", id, "check_several", params_str(dmz = c1$dmz, ppm = c1$ppm, scenario = c1$name),
        paste0(id, ".csv"))
}
# ppm=TRUE + overlap detected -> hits the trailing-comma bug in check_several.R
# and errors. Captured as the error message itself, not a normal output file,
# so task #8 can decide how to handle this known quirk with the reference in hand.
err <- tryCatch({ check_several(pat_overlap, dmz = 50, ppm = TRUE); NA_character_ },
    error = function(e) conditionMessage(e))
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(data.frame(error_message = err), file.path(dir_cs, paste0(id, ".csv")))
add_case("check_several", id, "check_several", params_str(dmz = 50, ppm = TRUE, scenario = "ppmTRUE_overlap"),
    paste0(id, ".csv"), notes = "known bug: trailing comma in paste0() call errors when ppm=TRUE and an overlap is found")

flush_manifest("check_several", dir_cs)
cat(sprintf("check_several: %d cases\n", n))

## =============================================================================
## 7. check_chemform
## =============================================================================
dir_cc <- file.path(OUT, "check_chemform"); dir.create(dir_cc, recursive = TRUE)
n <- 0

# 7a. broad pass, all 107 chemforms, default args
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(check_chemform(isotopes, chemforms), file.path(dir_cc, paste0(id, ".csv")))
add_case("check_chemform", id, "check_chemform", params_str(get_sorted = FALSE, get_list = FALSE, n_formulas = length(chemforms)),
    paste0(id, ".csv"), notes = "broad_pass_all_chemforms")

# 7b. edge cases: invalid formulas, nested brackets, isotope labels, D, empty string
edge_forms <- c(
    "C900Cl4H49", "O82394", "C8G500Zn9", "Br1", "6DBr9889", "",
    "C1(H2(O1)2)3", "D2H2O1", "[13]C2C8D2H8[15]N2N2[18]O2O3"
)
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(check_chemform(isotopes, edge_forms), file.path(dir_cc, paste0(id, ".csv")))
add_case("check_chemform", id, "check_chemform", params_str(get_sorted = FALSE, get_list = FALSE, n_formulas = length(edge_forms)),
    paste0(id, ".csv"), notes = "edge_cases: invalid formulas / nested brackets / isotope labels / D / empty string")

# 7c. get_sorted=TRUE, formulas with elements out of isotopes-table order
sort_forms <- c("H4C2O1", "O1H2C1", "N1H3C1O2")
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(check_chemform(isotopes, sort_forms, get_sorted = TRUE), file.path(dir_cc, paste0(id, ".csv")))
add_case("check_chemform", id, "check_chemform", params_str(get_sorted = TRUE, get_list = FALSE, n_formulas = length(sort_forms)),
    paste0(id, ".csv"), notes = "get_sorted=TRUE reorders elements to isotopes-table order")

# 7d. get_list=TRUE, flattened to long format (formula, element, count)
listed <- check_chemform(isotopes, pool10, get_list = TRUE)
long_rows <- lapply(names(listed), function(nm) {
    v <- listed[[nm]]
    if (length(v) == 0) return(data.frame(formula = nm, element = NA_character_, count = NA_integer_))
    data.frame(formula = nm, element = names(v), count = as.integer(v))
})
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(do.call(rbind, long_rows), file.path(dir_cc, paste0(id, ".csv")))
add_case("check_chemform", id, "check_chemform", params_str(get_sorted = FALSE, get_list = TRUE, n_formulas = length(pool10)),
    paste0(id, ".csv"), notes = "get_list=TRUE flattened to long format (formula,element,count)")

flush_manifest("check_chemform", dir_cc)
cat(sprintf("check_chemform: %d cases\n", n))

## =============================================================================
## 8. Formula-algebra helpers: mergeform, subform, multiform, check_ded
## =============================================================================
dir_fa <- file.path(OUT, "formula_algebra"); dir.create(dir_fa, recursive = TRUE)
n <- 0

# real adduct formula fragments, for realistic algebra inputs
add_frags <- unique(as.character(adducts$Formula_add[adducts$Formula_add != "FALSE" & nchar(as.character(adducts$Formula_add)) > 0]))
ded_frags <- unique(as.character(adducts$Formula_ded[adducts$Formula_ded != "FALSE" & nchar(as.character(adducts$Formula_ded)) > 0]))

# mergeform: formula2 must be scalar (see @param docs) - formula1 vectorizes.
# Loop over a few real adduct addition fragments, one .Call-level call each.
f1_pool <- pool10[1:5]
merge_rows <- lapply(head(add_frags, 3), function(frag) {
    data.frame(formula1 = f1_pool, formula2 = frag, result = mergeform(f1_pool, frag), stringsAsFactors = FALSE)
})
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(do.call(rbind, merge_rows), file.path(dir_fa, paste0(id, ".csv")))
add_case("formula_algebra", id, "mergeform", params_str(n_pairs = length(f1_pool) * 3), paste0(id, ".csv"))

# mergeform: D-label and isotope-label specific pairs (formula2 scalar per call)
f1m2 <- c("C10[13]C2H10Cl10", "D2C1H2O1", "C1H4")
f2m2 <- c("C2H5Na1", "[13]C1H2", "D4")
merge_rows2 <- lapply(seq_along(f2m2), function(i) {
    data.frame(formula1 = f1m2[i], formula2 = f2m2[i], result = mergeform(f1m2[i], f2m2[i]), stringsAsFactors = FALSE)
})
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(do.call(rbind, merge_rows2), file.path(dir_fa, paste0(id, ".csv")))
add_case("formula_algebra", id, "mergeform", params_str(n_pairs = length(f2m2), scenario = "D_and_isotope_labels"), paste0(id, ".csv"))

# subform: formula2/deduct is scalar too (same single-scalar parsing pattern
# as mergeform's formula2). Loop per scalar fragment; formula1 vectorizes.
sub_rows <- lapply(head(ded_frags, 3), function(frag) {
    ok <- !check_ded(pool10[1:5], frag)  # only keep pairs where subtraction is valid
    f1 <- pool10[1:5][ok]
    if (length(f1) == 0) return(NULL)
    data.frame(formula1 = f1, formula2 = frag, result = subform(f1, frag), stringsAsFactors = FALSE)
})
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(do.call(rbind, sub_rows), file.path(dir_fa, paste0(id, ".csv")))
add_case("formula_algebra", id, "subform", params_str(n_fragments = 3), paste0(id, ".csv"), notes = "valid_subtractions")

# subform: the two failure-message paths (element missing / count insufficient),
# one scalar formula2 per call
n <- n + 1; id <- sprintf("case_%03d", n)
fail_rows <- rbind(
    data.frame(formula1 = "C10H10", formula2 = "Cl1", result = subform("C10H10", "Cl1"), stringsAsFactors = FALSE),
    data.frame(formula1 = "C10H10", formula2 = "C20H10", result = subform("C10H10", "C20H10"), stringsAsFactors = FALSE)
)
write_precise(fail_rows, file.path(dir_fa, paste0(id, ".csv")))
add_case("formula_algebra", id, "subform", params_str(n_pairs = 2, scenario = "failure_message_paths"),
    paste0(id, ".csv"), notes = "returns a human-readable failure string instead of erroring")

# multiform: factor sweep on a curated formula
f_multi <- rep("C10[13]C2H10Cl10", 4)
facts <- c(1, 2, 3, 0.5)
res <- multiform(f_multi, facts)
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(data.frame(formula_in = f_multi, fact = facts, result = res), file.path(dir_fa, paste0(id, ".csv")))
add_case("formula_algebra", id, "multiform", params_str(n_pairs = length(res)), paste0(id, ".csv"))

# check_ded: deduct is scalar (same pattern); formulas vectorizes.
# Includes D-label normalisation (gsub D->H here, NOT [2]H like the other
# helpers - deliberate asymmetry to preserve in the port).
f_ded <- c("C8H4Cl2", "C10H16O2", "C3H10", "C4D10")
res <- check_ded(f_ded, "C4H10")
n <- n + 1; id <- sprintf("case_%03d", n)
write_precise(data.frame(formula = f_ded, deduct = "C4H10", warning = res), file.path(dir_fa, paste0(id, ".csv")))
add_case("formula_algebra", id, "check_ded", params_str(n_formulas = length(f_ded)), paste0(id, ".csv"),
    notes = "D normalised to plain H here (not [2]H like the other helpers) - deliberate asymmetry to preserve")

flush_manifest("formula_algebra", dir_fa)
cat(sprintf("formula_algebra: %d cases\n", n))

## =============================================================================
## Metadata: environment fingerprint, needed so #5's C extension can match
## compiler flags and #4/#9 can flag drift from a different R/gcc version.
## =============================================================================
gcc_version <- tryCatch(system("gcc --version", intern = TRUE)[1], error = function(e) NA_character_)
cflags <- tryCatch(system("R CMD config CFLAGS", intern = TRUE), error = function(e) NA_character_)
meta <- data.frame(
    key = c("r_version", "envipat_version", "gcc_version", "cflags", "platform"),
    value = c(R.version.string, as.character(packageVersion("enviPat")), gcc_version, cflags, R.version$platform),
    stringsAsFactors = FALSE
)
write.csv(meta, file.path(OUT, "metadata.csv"), row.names = FALSE)

cat("== done ==\n")
