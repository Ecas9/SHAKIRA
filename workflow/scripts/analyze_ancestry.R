#!/usr/bin/env Rscript
###############################################################################
# analyze_ancestry.R
#
# Post-process the SHAKIRA supervised ADMIXTURE run:
#   1. read prunedData.{Q,fam,pop}; split into 1KG reference vs study cell lines
#   2. map each Q column -> super-population EMPIRICALLY (from the labeled 1KG
#      individuals), so we never rely on ADMIXTURE's column order
#   3. plot the study cell lines' ancestry proportions (stacked bar)
#   4. join to the CCLE ancestry master reference (by normalized name / match_keys)
#      and score concordance: dominant-super-pop agreement + per-super-pop r / RMSE
#
# Outputs (under OUT_DIR):
#   ancestry_study_proportions.csv   per cell line: AFR/EUR/EAS/SAS/AMR + dominant
#   ancestry_concordance.csv         per matched cell line: ours vs reference
#   concordance_summary.txt          headline metrics
#   plot_ancestry_stacked.pdf        study ancestry stacked bars
#   plot_concordance_scatter.pdf     ours vs reference %, faceted by super-pop
#   plot_dominant_confusion.pdf      dominant super-pop confusion matrix
#
# Run:  Rscript analyze_ancestry.R         # uses the paths in the CONFIG block
#   or: Rscript analyze_ancestry.R <run_dir> <ref_csv> <out_dir>
###############################################################################

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(stringr); library(ggplot2); library(tibble)
})

# ============================== CONFIG =======================================
args <- commandArgs(trailingOnly = TRUE)
RUN_DIR <- if (length(args) >= 1) args[[1]] else "results"                       # normally supplied by Snakemake
REF_CSV <- if (length(args) >= 2) args[[2]] else "Validation/ccle_ancestry_master.csv"
OUT_DIR <- if (length(args) >= 3) args[[3]] else file.path(RUN_DIR, "analysis")

K       <- if (length(args) >= 4) as.integer(args[[4]]) else 5   # supervised K (optional 4th arg)
# SUPERPOPS is derived below from the reference labels actually present in the .pop (K-agnostic:
# AFR/AMR/EUR for a K=3 run, all five for K=5) — never a hardcoded list, so a K=3 panel run does
# not try to map empty EAS/SAS columns and crash.
# Files produced by the pipeline (ADMIXTURE writes <prefix>.K.Q next to where it ran):
Q_FILE   <- file.path(RUN_DIR, "admixture", sprintf("prunedData.%d.Q", K))
FAM_FILE <- file.path(RUN_DIR, "pruned", "prunedData.fam")
POP_FILE <- file.path(RUN_DIR, "pruned", "prunedData.pop")

SUPERPOP_COLORS <- c(AFR = "#E69F00", AMR = "#56B4E9", EAS = "#009E73",
                     EUR = "#D55E00", SAS = "#CC79A7")
# =============================================================================

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
norm_key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))   # JHH-1 -> JHH1
stopifnot(file.exists(Q_FILE), file.exists(FAM_FILE), file.exists(POP_FILE), file.exists(REF_CSV))

# ---- 1. load Q / fam / pop (all row-aligned, in .fam order) -----------------
Q   <- as.matrix(read.table(Q_FILE))
fam <- read.table(FAM_FILE, stringsAsFactors = FALSE)
iid <- fam[[2]]                                   # --double-id => FID==IID==sample
pop <- trimws(readLines(POP_FILE))
if (!(nrow(Q) == length(iid) && length(iid) == length(pop)))
  stop(sprintf("row mismatch: Q=%d fam=%d pop=%d", nrow(Q), length(iid), length(pop)))
if (ncol(Q) != K) stop(sprintf("Q has %d columns but K=%d", ncol(Q), K))

is_ref   <- nzchar(pop)                           # labeled 1KG = reference; blank = study cell line
cat(sprintf("[info] %d individuals: %d reference (labeled), %d study (blank)\n",
            length(iid), sum(is_ref), sum(!is_ref)))

# super-pops actually present in the reference labels (K-agnostic). Must equal K and Q's columns.
SUPERPOPS <- sort(unique(pop[is_ref]))
if (length(SUPERPOPS) != K)
  stop(sprintf("found %d distinct reference super-pop(s) [%s] but K=%d — pass the matching K",
               length(SUPERPOPS), paste(SUPERPOPS, collapse = ","), K))
cat("[info] reference super-pops:", paste(SUPERPOPS, collapse = ", "), "\n")

# ---- 2. EMPIRICAL Q-column -> super-pop map (from labeled 1KG individuals) ---
# For each super-pop, its reference individuals sit at ~1.0 in exactly one column;
# that column IS that super-pop. Derived from data, so column order is irrelevant.
ref_lab <- pop[is_ref]
colmeans_by_pop <- sapply(SUPERPOPS, function(s) colMeans(Q[is_ref, , drop = FALSE][ref_lab == s, , drop = FALSE]))
# colmeans_by_pop: K rows (columns) x length(SUPERPOPS); pick the best column per super-pop
col_for_pop <- apply(colmeans_by_pop, 2, which.max)
if (length(unique(col_for_pop)) != K)
  stop("Q-column -> super-pop mapping is not 1:1; inspect colmeans_by_pop:\n",
       paste(capture.output(print(round(colmeans_by_pop, 3))), collapse = "\n"))
col_label <- character(K); col_label[col_for_pop] <- names(col_for_pop)
cat("[info] Q column -> super-pop map: ", paste0("Q", seq_len(K), "=", col_label, collapse = "  "), "\n")
colnames(Q) <- col_label
Q <- Q[, SUPERPOPS]                               # reorder columns to AFR,AMR,EAS,EUR,SAS

# ---- 3. study cell-line ancestry table (percentages) ------------------------
studyP <- Q[!is_ref, , drop = FALSE] * 100          # matrix, cols = SUPERPOPS, in %
study <- tibble(cell_line = iid[!is_ref]) %>%
  bind_cols(as_tibble(studyP)) %>%
  mutate(dominant     = SUPERPOPS[max.col(studyP, ties.method = "first")],
         dominant_pct = apply(studyP, 1, max),
         normkey      = norm_key(cell_line)) %>%
  arrange(dominant, desc(dominant_pct))
write_csv(study %>% select(-normkey), file.path(OUT_DIR, "ancestry_study_proportions.csv"))

# stacked bar of study ancestry
study_long <- study %>%
  mutate(cell_line = factor(cell_line, levels = cell_line)) %>%
  pivot_longer(all_of(SUPERPOPS), names_to = "superpop", values_to = "pct")
gbar <- ggplot(study_long, aes(cell_line, pct, fill = superpop)) +
  geom_col(width = 1) +
  scale_fill_manual(values = SUPERPOP_COLORS) +
  coord_flip() + labs(x = NULL, y = "Ancestry (%)", fill = "Super-pop",
                      title = sprintf("Supervised ADMIXTURE ancestry — %d cell lines", nrow(study))) +
  theme_minimal(base_size = 7) + theme(legend.position = "top")
ggsave(file.path(OUT_DIR, "plot_ancestry_stacked.pdf"), gbar,
       width = 7, height = max(4, nrow(study) * 0.13), limitsize = FALSE)

# ---- 4. join to the CCLE master reference by normalized key -----------------
ref <- read_csv(REF_CSV, show_col_types = FALSE)
req <- c("cLine_ID", "dominant_superpop", "match_keys", SUPERPOPS)
miss_cols <- setdiff(req, names(ref))
if (length(miss_cols)) stop("ccle_ref_csv missing required column(s): ", paste(miss_cols, collapse = ", "))
# expand every match key (and the normalized cLine_ID) into a long lookup
ref_keys <- ref %>%
  mutate(.row = row_number(),
         keys = paste(norm_key(cLine_ID),
                      toupper(ifelse(is.na(match_keys), "", match_keys)), sep = "|")) %>%
  separate_rows(keys, sep = "\\|") %>%
  filter(keys != "") %>%
  distinct(keys, .keep_all = TRUE) %>%
  dplyr::select(.row, keys)

ref_ref <- ref %>% transmute(.row = row_number(),
                             ref_dominant = dominant_superpop,
                             ref_ethnicity = if ("Ethnicity" %in% names(ref)) Ethnicity else NA_character_)
# attach the reference proportion columns for the present super-pops as ref_<SUPERPOP>
ref_ref <- bind_cols(ref_ref, setNames(ref[SUPERPOPS], paste0("ref_", SUPERPOPS)))

matched <- study %>%
  left_join(ref_keys, by = c("normkey" = "keys")) %>%
  left_join(ref_ref, by = ".row")

n_match <- sum(!is.na(matched$.row))
unmatched <- matched %>% filter(is.na(.row)) %>% pull(cell_line)
cat(sprintf("[info] matched %d / %d study cell lines to the reference\n", n_match, nrow(study)))
if (length(unmatched)) cat("[info] unmatched:", paste(unmatched, collapse = ", "), "\n")

conc <- matched %>% filter(!is.na(.row)) %>%
  dplyr::select(cell_line, dominant, ref_dominant, ref_ethnicity,
         all_of(SUPERPOPS), all_of(paste0("ref_", SUPERPOPS))) %>%
  mutate(dominant_match = dominant == ref_dominant)
write_csv(conc, file.path(OUT_DIR, "ancestry_concordance.csv"))

# ---- 5. concordance metrics -------------------------------------------------
dom_agree <- mean(conc$dominant_match)
permetric <- lapply(SUPERPOPS, function(s) {
  o <- conc[[s]]; r <- conc[[paste0("ref_", s)]]
  tibble(superpop = s, pearson_r = suppressWarnings(cor(o, r)),
         rmse = sqrt(mean((o - r)^2)), mae = mean(abs(o - r)))
}) %>% bind_rows()

summ <- c(
  sprintf("Matched cell lines:            %d / %d", n_match, nrow(study)),
  sprintf("Dominant super-pop agreement:  %.1f%% (%d/%d)",
          100 * dom_agree, sum(conc$dominant_match), nrow(conc)),
  sprintf("Mean abs error (all props):    %.2f%%", mean(permetric$mae)),
  "", "Per super-population:",
  capture.output(print(as.data.frame(permetric), row.names = FALSE, digits = 3)))
writeLines(summ, file.path(OUT_DIR, "concordance_summary.txt"))
cat("\n", paste(summ, collapse = "\n"), "\n", sep = "")

# ---- 6. concordance plots ---------------------------------------------------
sc_long <- conc %>%
  dplyr::select(cell_line, all_of(SUPERPOPS), all_of(paste0("ref_", SUPERPOPS))) %>%
  pivot_longer(-cell_line, names_to = "k", values_to = "v") %>%
  mutate(which = ifelse(str_starts(k, "ref_"), "ref", "ours"),
         superpop = str_remove(k, "ref_")) %>%
  dplyr::select(cell_line, superpop, which, v) %>%
  pivot_wider(names_from = which, values_from = v)
gsc <- ggplot(sc_long, aes(ref, ours, color = superpop)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(alpha = 0.6, size = 1.2) +
  facet_wrap(~superpop) + scale_color_manual(values = SUPERPOP_COLORS, guide = "none") +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(x = "Reference (%)", y = "SHAKIRA Estimates (%)",
       title = "Per-super-population concordance vs CCLE reference") +
  theme_bw(base_size = 20)
ggsave(file.path(OUT_DIR, "plot_concordance_scatter.pdf"), gsc, width = 7, height = 5)

cm1 <- conc %>% dplyr::count(ref_dominant, dominant) %>%
  complete(ref_dominant = SUPERPOPS, dominant = SUPERPOPS, fill = list(n = 0))
gcm <- ggplot(cm1, aes(ref_dominant, dominant, fill = n)) +
  geom_tile(color = "white") + geom_text(aes(label = n), size = 3) +
  scale_fill_gradient(low = "white", high = "#3182bd") +
  labs(x = "Reference dominant super-pop", y = "SHAKIRA dominant super-pop",
       title = sprintf("Dominant super-pop \n confusion table (%.0f%% agree)", 100 * dom_agree)) +
  theme_minimal(base_size = 15) + theme(panel.grid = element_blank(),plot.title = element_text(hjust = 0.5))
ggsave(file.path(OUT_DIR, "plot_dominant_confusion.pdf"), gcm, width = 5.5, height = 5)

cat("\n[done] wrote results + plots to:", OUT_DIR, "\n")
