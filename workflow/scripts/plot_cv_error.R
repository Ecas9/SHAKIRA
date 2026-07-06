#!/usr/bin/env Rscript
###############################################################################
# plot_cv_error.R
#
# Parse the ADMIXTURE cross-validation logs (the SHAKIRA admixture_cv rule) and plot
# CV error vs K to pick the best K (the minimum). Run AFTER the CV sweep.
#
# Inputs:  <run_dir>/admixture/cv/log.<K>.out  (one per K; ADMIXTURE prints one
#          line like "CV error (K=3): 0.45297"). K is parsed from the FILENAME,
#          not the line text — the legacy lambda-sweep logs all say K=3 internally.
# Outputs (under <out_dir>, default <run_dir>/analysis):
#          cv_error.csv         K,error (sorted by K)
#          plot_cv_error.pdf    error vs K, minimum marked in red
#
# Run:  Rscript plot_cv_error.R <run_dir> [out_dir]
#
# Kept separate from analyze_ancestry.R on purpose: the CV plot only needs the CV
# logs, so it still runs for a K=3 panel run even when the concordance analysis
# (which expects the CCLE master + all super-pops) does not apply.
###############################################################################
suppressPackageStartupMessages(library(ggplot2))

args    <- commandArgs(trailingOnly = TRUE)
RUN_DIR <- if (length(args) >= 1) args[[1]] else "results"   # normally supplied by Snakemake
OUT_DIR <- if (length(args) >= 2) args[[2]] else file.path(RUN_DIR, "analysis")
CV_DIR  <- file.path(RUN_DIR, "admixture", "cv")

files <- list.files(CV_DIR, pattern = "^log\\.[0-9]+\\.out$", full.names = TRUE)
if (length(files) == 0) {
  cat("[plot_cv_error] no log.<K>.out files in", CV_DIR, "- nothing to plot\n")
  quit(status = 0)
}

rows <- lapply(files, function(f) {
  K  <- as.integer(sub("^log\\.([0-9]+)\\.out$", "\\1", basename(f)))   # K from FILENAME
  ln <- grep("CV error", readLines(f, warn = FALSE), value = TRUE)
  if (length(ln) == 0) return(NULL)
  err <- as.numeric(sub(".*:[[:space:]]*([0-9.]+).*", "\\1", ln[[1]]))
  data.frame(K = K, error = err)
})
cv <- do.call(rbind, rows)
if (!is.null(cv)) cv <- cv[!is.na(cv$error), , drop = FALSE]   # drop logs whose CV value didn't parse
if (is.null(cv) || nrow(cv) == 0) {
  cat("[plot_cv_error] no parseable 'CV error' values found under", CV_DIR, "\n")
  quit(status = 0)
}
cv <- cv[order(cv$K), ]
best_K <- cv$K[which.min(cv$error)]

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(cv, file.path(OUT_DIR, "cv_error.csv"), row.names = FALSE)

g <- ggplot(cv, aes(K, error)) +
  geom_line() + geom_point(size = 2) +
  geom_point(data = cv[cv$K == best_K, , drop = FALSE], color = "red", size = 3.5) +
  scale_x_continuous(breaks = cv$K) +
  labs(x = "K", y = "Cross-validation error",
       title = sprintf("ADMIXTURE CV error (minimum at K=%d)", best_K)) +
  theme_bw(base_size = 11)
ggsave(file.path(OUT_DIR, "plot_cv_error.pdf"), g, width = 6, height = 4)

cat(sprintf("[plot_cv_error] %d K value(s); best K=%d (CV error %.5f)\n",
            nrow(cv), best_K, min(cv$error)))
cat("[plot_cv_error] wrote", file.path(OUT_DIR, "cv_error.csv"), "and plot_cv_error.pdf\n")
