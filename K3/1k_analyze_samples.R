#!/usr/bin/env Rscript
# Identify 1KG samples that are NOT >= `thr` (default 0.98) in ANY single ancestry at
# the chosen K -- i.e. the admixed ones -- and write them as the exclude list
# (FAMID<TAB>sample) consumed by exclude_admixed_1K_samples.sh (plink --remove) to build
# the high-quality / "pure-ancestry" 1KG reference panel.
#
# K-agnostic: instead of naming the K=3 columns (EUR/AMR/AFR), keep any row whose MAX
# ancestry proportion is < thr -- the same rule for any K.
#
# Usage: Rscript 1k_analyze_samples.R [K] [prefix] [base_dir] [des_dir] [thr]
#   K        supervised K used in the ADMIXTURE run   (default 3)
#   prefix   bed/.fam/.pop/.Q basename                (default Final_high_quality_default)
#   base_dir dir holding <prefix>.{fam,pop,<K>.Q}     (default /data/.../1K)
#   des_dir  dir to write 1K_98_co_exclude.txt        (default /data/.../K3_15_mill_filter)
#   thr      "pure ancestry" threshold                (default 0.98)
library(tidyverse)

args     <- commandArgs(trailingOnly = TRUE)
K        <- if (length(args) >= 1) as.integer(args[[1]]) else 3L
prefix   <- if (length(args) >= 2) args[[2]] else "Final_high_quality_default"
base_dir <- if (length(args) >= 3) args[[3]] else "/path/to/data/1K"
des_dir  <- if (length(args) >= 4) args[[4]] else "/path/to/data/K3_15_mill_filter"
thr      <- if (length(args) >= 5) as.numeric(args[[5]]) else 0.98

fam.df <- read_delim(file.path(base_dir, paste0(prefix, ".fam")), col_names = FALSE) %>%
    rename(sample = X2) %>%
    mutate(pop = read.table(file.path(base_dir, paste0(prefix, ".pop")))$V1)

Q.df <- read_delim(file.path(base_dir, paste0(prefix, ".", K, ".Q")), col_names = FALSE)

anc.df <- Q.df %>%
    mutate(sample  = fam.df$sample,
           pop     = fam.df$pop,
           max_anc = do.call(pmax, dplyr::select(Q.df, dplyr::all_of(paste0("X", seq_len(K)))))) %>%
    filter(max_anc < thr) %>%       # no single ancestry reaches thr -> admixed -> exclude
    mutate(FAMID = 0)

write_tsv(dplyr::select(anc.df, FAMID, sample),
          file.path(des_dir, "1K_98_co_exclude.txt"), col_names = FALSE)

message(sprintf("K=%d thr=%.2f: %d admixed sample(s) excluded -> %s/1K_98_co_exclude.txt",
                K, thr, nrow(anc.df), des_dir))
