# K3 high-quality-reference scripts (legacy — provenance only)

> **Not the maintained pipeline.** These are the original SLURM scripts that build the optional
> "<98%-ancestry" 1000 Genomes exclusion list (`1K_98_co_exclude.txt`) used by the K=3 run. They
> are kept for provenance and are **not portable**: they hardcode UAB Cheaha modules, conda env
> names, absolute `/scratch` + `/data/project` paths, and a personal `--mail-user`. Do not run
> them as-is. The maintained workflow is the Snakemake pipeline at the repo root (see the
> top-level `README.md`); it drops admixed subpopulations via `exclude_subpops` + per-superpop MAF
> and exposes the optional exclusion list through `reference.kg_exclude` in `config.yaml`.

Contents: `run_1k_highqual_cheaha.sh` -> `1k_analyze_samples.R` (build the exclusion list),
`1kGenome_filter.sh`, `exclude_admixed_1K_samples.sh`, `1k_filter.py`, `pop_jhu.py`.
