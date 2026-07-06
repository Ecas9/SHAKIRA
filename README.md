# SHAKIRA

**Superpopulation-level Haplotype Ancestry from K-population Inference of RNA Assays** — a Snakemake
workflow that estimates continental genetic ancestry for a set of samples by supervised ADMIXTURE
against a 1000 Genomes reference panel, then scores concordance against a known-ancestry reference
(CCLE). Study genotypes come from an upstream nf-core variant-calling run (`rnavar` for RNA-seq, the
default; `sarek` for DNA germline).

The pipeline was developed on UAB's Cheaha HPC and has been generalized so it can be run on any
SLURM cluster (or a single machine) by editing configuration only. Nothing site-specific is baked
into the workflow logic.

## What it does

Per-sample caller VCFs → normalize/QC → harmonize against the 1000G panel → merge → variant QC and
LD pruning → **supervised ADMIXTURE** (K reference superpopulations) → ancestry proportions +
concordance report + a cross-validation plot to choose K.

Default run: K=3 supervised with AFR/EUR/AMR reference superpopulations. The reference is restricted
to the superpopulations of interest, admixed subpopulations are dropped, and per-superpopulation MAF
is applied so ancestry-informative markers common in one group are retained (see
`workflow/scripts/ref_maf_per_superpop.sh`).

Each Snakemake rule is one pipeline stage; the per-sample ingest stage fans out into one job per
sample.

## Requirements

- **conda / mamba** (Miniforge recommended). All tools are pinned in `environment.yml` and the
  per-rule `envs/*.yaml`; Snakemake materializes them with `--use-conda`. No environment modules
  required.
- A **SLURM** cluster for the default profile, or use the bundled `profiles/local` to run on one
  machine.
- **Singularity/Apptainer + Nextflow** only for the optional upstream variant calling (nf-core
  rnavar/sarek). These run in their own containers; see `Validation/`.

Exact versions are in `TOOLS_VERSIONS.txt`; `refs/capture_versions.sh` records what is actually
installed in an env.

## Install

```bash
git clone <this-repo> && cd <this-repo>
conda env create -f environment.yml     # creates env "shakira"
conda activate shakira
```

The `shakira` env holds Snakemake and the SLURM executor plugin (the controller). The
per-rule scientific tools are installed on demand by `--use-conda`.

## Configure

All paths are built from three **anchors** in `config.yaml` under `path_anchors`. Set them once,
either by editing `config.yaml` or by exporting the matching environment variables (anchor and path
values also expand `${VARS}` and `~`):

| anchor        | env override     | what it points to                                              |
|---------------|------------------|----------------------------------------------------------------|
| `work_root`   | `ANCESTRY_WORK`  | writable area: caller output (study VCFs), all outputs, caller ref (FASTA/GTF) |
| `panel_dir`   | `ANCESTRY_PANEL` | directory holding the 1000G panel VCF                           |
| `gatk_bundle` | `ANCESTRY_GATK`  | GATK known-sites bundle (upstream variant calling only)        |

```bash
export ANCESTRY_WORK=/scratch/$USER/ancestry_work
export ANCESTRY_PANEL=/data/refs/1kg
export ANCESTRY_GATK=/data/refs/gatk_bundle/GRCh38
```

The shipped anchors are `/path/to/...` placeholders. A run whose anchors are still unset fails
immediately with a message naming the anchor to set, so it can never silently use the wrong data.
Analysis parameters (K, MAF, LD-pruning, superpopulations, etc.) are locked in the `params:` block
of `config.yaml`; edit there.

## Reference data

Provide and validate references once:

1. **1000G GRCh38 panel** `1kgenome.vcf.gz` + index (`.tbi`/`.csi`) under `panel_dir`. Bare contigs
   (1..22), bgzipped + indexed. Provenance is recorded in `config.yaml`.
2. **Ensembl FASTA + GTF v112** under `work_root/ref/...` (upstream calling only).
3. Run `bash refs/prep_refs.sh` (uses the same anchors). It hard-fails if the panel is missing, and
   renames the GATK known-sites from chr-prefixed to bare contigs, reheadered to the FASTA.

The 1000G sample→population table (`resources/1K_pops.txt`) and the long-range-LD exclusion regions
(`resources/high-LD-regions-hg38-GRCh38.txt`) are committed with the repo.

## Upstream variant calling (produce the study VCFs)

If you do not already have per-sample caller VCFs, run the nf-core pipeline (`Validation/`):

```bash
bash init_nfcore_singularity.sh                 # one-time: pull rnavar/sarek + stage containers
mkdir -p Validation/slurm_logs
sbatch Validation/run_rnavar.slurm              # RNA-seq (default); or run_sarek.slurm for DNA
```

`nfcore_hpc.config` carries the SLURM+Singularity resource settings (partition via `NFCORE_QUEUE`;
resource caps tuned to Cheaha, edit for your nodes). The runners are configured with the same
anchors. Set `caller:` in `config.yaml` to match (`rnavar` or `sarek`).

## Run

```bash
# SLURM (submits one job per rule via profiles/slurm)
snakemake -s workflow/Snakefile --use-conda --profile profiles/slurm
# or the convenience wrapper (activates the conda env first):
bash run_snakemake.sh
bash run_snakemake.sh -n            # dry run: preview the DAG
bash run_snakemake.sh --unlock      # release a stale lock

# Single machine, no scheduler:
SMK_PROFILE=profiles/local bash run_snakemake.sh
```

## Outputs (under `work_root`/`outdir`, default `<work_root>/amr`)

- `admixture/prunedData.<K>.Q` — supervised ancestry proportions (rows aligned to `prunedData.fam`)
- `analysis/ancestry_study_proportions.csv` — per-sample AFR/EUR/AMR(/…) proportions + dominant call
- `analysis/concordance_summary.txt`, `ancestry_concordance.csv` — agreement vs the CCLE reference
- `analysis/plot_cv_error.pdf` — cross-validation error vs K (choose K at the minimum)
- `analysis/plot_ancestry_stacked.pdf`, `plot_concordance_scatter.pdf`, `plot_dominant_confusion.pdf`

## Porting to another cluster

The workflow logic is site-agnostic; adapt only configuration:

- **Paths**: set the three anchors (above). Nothing else contains a path to edit.
- **Partitions / account**: set them in `profiles/slurm/config.yaml` only. By default no
  `--partition`/`--account` is sent, so jobs use your cluster's default partition. Uncomment
  `default-resources.slurm_partition` / `slurm_account` to pin them, and the optional `set-resources`
  block to route the three short rules to a fast/express queue. The Snakefile itself carries no
  partitions.
- **Walltime / memory**: each rule declares portable `runtime` (minutes) and `mem_mb`; the SLURM
  executor maps these to `--time`/`--mem`. Adjust in the rule or via profile `set-resources` if your
  nodes are smaller.
- **No SLURM**: use `profiles/local` (set `cores:` to the machine).
- **No Lmod**: the `module load` calls in the wrappers are guarded and skipped automatically; just
  have conda (and, for upstream calling, Nextflow + Singularity) on `PATH`.
- **nf-core resources**: edit `Validation/nfcore_hpc.config` `resourceLimits` and the withLabel/
  withName sizes to your hardware; set the partition with `NFCORE_QUEUE`.

## Repository layout

```
config.yaml               single source of truth (anchors, params, caller toggle, env map)
environment.yml           controller conda env (snakemake + slurm plugin)
run_snakemake.sh          launcher (CONDA_ENV / SMK_PROFILE env knobs)
init_nfcore_singularity.sh   one-time nf-core bootstrap
envs/                     per-rule conda envs
profiles/  slurm/  local/    execution profiles (site config lives here)
refs/      prep_refs.sh  capture_versions.sh
resources/                committed reference tables (1K_pops, high-LD regions)
workflow/  Snakefile  scripts/    the workflow + its helper scripts
Validation/               nf-core rnavar/sarek runners + resource config + sample-sheet builder
cheaha/  K3/              legacy scripts, kept for provenance only (see their README.md)
```

## Reproducibility

Tool versions are pinned in `environment.yml`, `envs/*.yaml`, and documented in `TOOLS_VERSIONS.txt`
(as-run vs standardized). nf-core pipeline revisions are pinned by `-r` (rnavar 1.3.0, sarek 3.9.0);
their internal tool versions are fixed by the pipeline release and its Singularity containers.

## Legacy scripts

`cheaha/` (the original single-job standalone driver) and `K3/` (the optional high-quality-reference
builder) are the pre-Snakemake scripts this workflow was ported from. They are kept for provenance,
are **not maintained**, and are **not portable** (they hardcode Cheaha modules and absolute paths).
Use the Snakemake workflow above. See `cheaha/README.md` and `K3/README.md`.

## Citation

If you use SHAKIRA, please cite:

Castillo E II, Herring B, Zaibaq F, Elhussin I, Concors S, Gillis A, Yates C, Rose JB, Guenter R.
*SHAKIRA: Superpopulation-level Haplotype Ancestry from K-population Inference of RNA Assays.*
(manuscript, 2026).

## Data provenance

- 1000 Genomes GRCh38 panel: 20190312 biallelic SNV+INDEL release (Lowy-Gallego et al. 2019,
  *Wellcome Open Research* 4:50).
- CCLE ancestry reference: Dutil et al. 2019 (used by `Validation/map_srr_to_ancestry.py` and the
  concordance step).
