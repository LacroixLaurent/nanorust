#!/usr/bin/env Rscript
# Compare a nanorust alldata.rds with the R pipeline's one.
# Usage: compare_alldata.R new.rds ref.rds [chrom1,chrom2,...]
# Rows are matched after ordering by (chrom, read_id, flag, start); doubles use a relative tolerance.
suppressMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
new <- readRDS(args[1])
ref <- readRDS(args[2])
if (length(args) >= 3) {
  chroms <- strsplit(args[3], ",")[[1]]
  ref <- ref %>% filter(chrom %in% chroms)
  new <- new %>% filter(chrom %in% chroms)
}
tol <- 1e-12
ok <- TRUE
report <- function(what, pass, detail = "") {
  cat(sprintf("%-38s %s %s\n", what, if (pass) "OK  " else "FAIL", detail))
  if (!pass) ok <<- FALSE
}

report("columns", identical(names(new), names(ref)), paste(names(new), collapse = ","))
report("column types", identical(sapply(new, typeof), sapply(ref, typeof)))
report("class", identical(class(new), class(ref)))
report("chrom levels", identical(levels(new$chrom), levels(ref$chrom)))
report("strand levels", identical(levels(new$strand), levels(ref$strand)))
report("nrow", nrow(new) == nrow(ref), sprintf("%d vs %d", nrow(new), nrow(ref)))

ord <- function(d) d %>% arrange(chrom, read_id, flag, start)
new <- ord(new); ref <- ord(ref)
key <- function(d) paste(d$read_id, d$flag, d$chrom, d$start)
only_new <- setdiff(key(new), key(ref)); only_ref <- setdiff(key(ref), key(new))
report("same mappings", length(only_new) == 0 && length(only_ref) == 0,
       sprintf("only_new=%d only_ref=%d", length(only_new), length(only_ref)))
if (length(only_new)) print(head(new[key(new) %in% only_new, 1:6]))
if (length(only_ref)) print(head(ref[key(ref) %in% only_ref, 1:6]))

common <- intersect(key(new), key(ref))
n <- new[match(common, key(new)), ]; r <- ref[match(common, key(ref)), ]
for (col in c("read_id", "flag", "chrom", "strand", "start", "end"))
  report(col, identical(n[[col]], r[[col]]))
reldiff <- function(a, b) { d <- abs(a - b) / pmax(abs(b), 1e-300); d[a == b] <- 0; max(c(0, d)) }
for (col in c("med_signal", "med_signalbin")) {
  rd <- reldiff(n[[col]], r[[col]])
  report(col, rd <= tol, sprintf("max rel diff %.3g, bit-identical %.4f", rd, mean(n[[col]] == r[[col]])))
}
nb_same <- mapply(function(a, b) nrow(a) == nrow(b) && identical(a$positions, b$positions), n$signalbin, r$signalbin)
report("signalbin positions", all(nb_same), sprintf("%d mismatching", sum(!nb_same)))
if (all(nb_same)) {
  sb <- mapply(function(a, b) reldiff(a$signalB, b$signalB), n$signalbin, r$signalbin)
  bit <- mapply(function(a, b) all(a$signalB == b$signalB), n$signalbin, r$signalbin)
  report("signalbin signalB", max(sb) <= tol, sprintf("max rel diff %.3g, bit-identical rows %.4f", max(sb), mean(bit)))
}
sa <- function(x) lapply(attributes(x), function(a) if (is.numeric(a)) class(a) else a)
report("nested tibble attributes", identical(sa(n$signalbin[[1]]), sa(r$signalbin[[1]])))
cat(if (ok) "ALL CHECKS PASSED\n" else "SOME CHECKS FAILED\n")
quit(status = if (ok) 0 else 1)
