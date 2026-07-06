# Legacy standalone runner (provenance — not the maintained pipeline)

> **This directory is kept for provenance only.** It holds the original pre-Snakemake
> driver the published `workflow/Snakefile` was ported from. It is **not maintained** and
> is **not portable** (it hardcodes UAB Cheaha modules, conda env names, and absolute
> `/scratch` + `/data/project` paths). To actually run the pipeline, use the Snakemake
> workflow at the repo root — see the top-level `README.md`.
>
> The active helper scripts and reference tables that used to live here have moved:
> `ingest_one.sh`, `ref_maf_per_superpop.sh`, `pop.py`, `analyze_ancestry.R`,
> `plot_cv_error.R` → `workflow/scripts/`; `1K_pops.txt`,
> `high-LD-regions-hg38-GRCh38.txt` → `resources/`. The copies referenced below no
> longer exist in this folder.

`run_ancestry_cheaha.sh` runs the whole RNAvar → ADMIXTURE ancestry pipeline as **one
SLURM job**, using Cheaha's **pre-existing modules + conda envs** (`admix`, `plink2`) —
no Snakemake, no per-rule conda-env files. It mirrors the logic of `../workflow/Snakefile`.

**Reference build: GRCh38 only.** The rnavar study VCFs and the 1000G panel are both
GRCh38 with bare contigs (1..22), so there is no liftover and no per-chromosome
fan-out — the flow is linear. (An hg19/liftover variant is available in the Snakemake
version via `ref_build: hg19`; it is intentionally not in this script.)

Files here:
- `run_ancestry_cheaha.sh` — the driver (edit the EDIT block, then `sbatch` it). This is the
  canonical driver (formerly `run_ancestry_cheaha_fix_MAF.sh`); older versions are in `archive/`.
- `pop.py` — builds the supervised `.pop` labels (copied from the repo)
- `ingest_one.sh`, `ref_maf_per_superpop.sh` — Stage 1 / Stage 4a helpers
- `analyze_ancestry.R`, `plot_cv_error.R` — post-ADMIXTURE analysis (see `run_analysis_cheaha.sh`)
- `1K_pops.txt`, `high-LD-regions-hg38-GRCh38.txt` — reference tables

**Snakemake migration (repo root):** `config.yaml` (paths + locked params + caller toggle),
`environment.yml` / `envs/*.yaml` (conda envs replacing the module loads), `refs/prep_refs.sh`
(reference validation + GATK rename), `init_nfcore_singularity.sh` (nextflow → nf-core → rnavar/sarek),
`TOOLS_VERSIONS.txt`, and `PREFLIGHT_snakemake.md` (status + remaining gaps).

## Prerequisites
- Modules on PATH: `BCFtools`, `PLINK/1.90-foss-2016a`, `Anaconda3`
- Your conda envs: `plink2` (provides `plink2`), `admix` (provides `admixture`)
- The GRCh38 panel VCF **bgzipped + tabix-indexed** (`.tbi`/`.csi`)

## Run
```bash
# 1) edit the SBATCH header (account/partition/mail) and the EDIT block, then:
sbatch cheaha/run_ancestry_cheaha.sh

# 2) watch it
squeue -u $USER
tail -f rnavar_admixture_*.out

# 3) results
#    supervised proportions:  $OUTDIR/admixture/prunedData.5.Q   (rows match prunedData.fam)
#    pick K from CV error:     grep -h CV $OUTDIR/admixture/cv/log.*.out
```

## Parameters to set (top of the script, "EDIT BLOCK")
| variable | what to set |
|---|---|
| `REPO_DIR` | your checkout path (locates `pop.py`, `1K_pops.txt`, pop table) |
| `INPUT_MODE` | `discover` = `find` all VCFs under `RNAVAR_DIR` and **write** `SAMPLE_SHEET`; `sheet` = **read** an existing `SAMPLE_SHEET` |
| `RNAVAR_DIR` / `SAMPLE_SHEET` | rnavar output dir to glob, and the 2-col `sample <tab> vcf` TSV (written in discover mode, read in sheet mode) |
| `KG_GRCH38` | combined GRCh38 1000G panel VCF (must be tabix-indexed) |
| `KG_PANEL` | 1KG sample→pop panel (`step_4/1K_pops.txt`). pop→super-pop mapping is built into `pop.py` — no external table needed |
| `OUTDIR` | where all outputs go |
| `K_SUPERVISED`, `K_MIN`, `K_MAX` | supervised K and the CV sweep range |
| `MAX_MISSING`, `LD_WINDOW/STEP/R2` | QC + LD-pruning thresholds |
| `CONDA_ENV_PLINK2`, `CONDA_ENV_ADMIX` | your env names if not `plink2` / `admix` |

Also edit `#SBATCH --mail-user` (and partition/time/cpus/mem) in the header.
