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

### Two settings worth checking before your first run

These are the two ways a run silently produces nothing, so they are called out separately. Both are
handled at ingest, in `workflow/scripts/ingest_one.sh`.

**`reference.contig_style`** - the contig naming used by *your 1000G panel*: `"nochr"` (`1`, `2`, ...
`22`, the default and what the 1000G releases use) or `"chr"` (`chr1`, `chr2`, ...). Your study VCFs
do **not** have to match: ingest detects each study VCF's own naming from its first data record and
converts it. Set this to whatever the panel uses. A study/panel contig mismatch produces an empty
intersection with no error, which is the single most common way to get a run that "works" and
returns nothing.

**`params.ingest_filters`** - the `bcftools --apply-filters` list applied to each study VCF. The
default `"PASS,."` keeps records marked `PASS` *and* records whose `FILTER` column is unset (`.`),
which is what GATK HaplotypeCaller writes when no `VariantFiltration` / `CNNScoreVariants` step was
run. A bare `"PASS"` yields **zero** variants for those inputs. Use `"PASS"` to require an explicit
pass, or e.g. `"PASS,.,LowQual"` to be more permissive.

Ingest also forces each single-sample VCF's internal sample name to the pipeline's sample ID, so the
cohort merge, the QC table and the ADMIXTURE output rows all agree even when the caller wrote a
run-level name (nf-core commonly writes e.g. `patient1_B30` rather than `B30`), and two inputs can
never collide on a shared internal name.

## Reference data

Provide and validate references once:

1. **1000G GRCh38 panel** `1kgenome.vcf.gz` + index (`.tbi`/`.csi`) under `panel_dir`. Bare contigs
   (1..22), bgzipped + indexed. Provenance is recorded in `config.yaml`.
2. **Ensembl FASTA + GTF v112** under `work_root/ref/...` (upstream calling only).
3. Run `bash refs/prep_refs.sh` (uses the same anchors). It hard-fails if the panel is missing, and
   renames the GATK known-sites from chr-prefixed to bare contigs, reheadered to the FASTA.

The 1000G sample→population table (`resources/1K_pops.txt`) and the long-range-LD exclusion regions
(`resources/high-LD-regions-hg38-GRCh38.txt`) are committed with the repo.

## Try it on example data

No access to a restricted cohort is needed to test SHAKIRA. `resources/get_example_data.sh`
downloads six **public CCLE RNA-seq runs** (SRA study SRP186687) whose DNA-based ancestry is already
published in `Validation/ccle_ancestry_master.csv`, so the correct answer is known in advance:

```bash
bash resources/get_example_data.sh --list          # show the six cell lines and expected ancestry
bash resources/get_example_data.sh                 # ~1.5 GB, a few minutes
bash resources/get_example_data.sh --full          # complete runs (~90 GB) if you want full depth
```

By default it uses HTTP range requests to pull only the first few hundred MB of each gzipped FASTQ
and truncates to whole records, giving ~2M read pairs per sample. That is shallow compared with a
full run but still yields tens of thousands of autosomal SNPs per sample - well above the
`params.qc.min_snps` gate - so the ancestry calls come out correct.

It writes `example_data/samplesheet.csv` (ready for nf-core/rnavar) and
`example_data/expected_ancestry.csv`. Run the two pipeline stages, then compare:

| sample | cell line | expected |
|---|---|---|
| Raji, P3HR1, PLCPRF5 | Burkitt lymphoma / hepatoma | >95% AFR |
| MDAMB361, WM2664, VMCUB1 | breast / melanoma / bladder | >95% EUR |

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

## Quality control

**Sample-level input gate.** After ingest and *before* the cohort merge, every sample is scored on
four metrics and samples that cannot support an ancestry call are excluded with a recorded reason.
A library with too few callable SNPs, or with coverage too low for confident genotypes, would
otherwise be carried all the way to ADMIXTURE and return an ancestry vector driven by a handful of
sites. Thresholds live in `config.yaml` under `params.qc`:

| key | default | what it catches |
|---|---|---|
| `min_snps` | 10000 | too few retained biallelic autosomal SNPs to place the sample |
| `min_median_dp` | 0 (off) | median `FORMAT/DP`; **reported, not gated** for RNA-seq - see below |
| `het_range` | `[0.05, 0.70]` | heterozygous fraction outside this = contamination (high) or clonal/LOH (low) |
| `min_ts_tv` | 1.50 | Ti/Tv collapse toward 0.5 indicates noise rather than real germline SNVs |
| `max_fail_frac` | 0.50 | aborts the run if most samples fail, which means an upstream problem |

Outputs: `qc/sample_qc.tsv` (every sample, every metric, PASS/FAIL + reason), `qc/samples_pass.txt`,
`qc/samples_fail.txt`. Set a threshold to `0` to disable that check.

**Why `min_median_dp` is off by default.** `FORMAT/DP` at a called RNA-seq SNP tracks transcript
abundance, not library quality: most called sites sit in moderately expressed transcripts, so the
per-sample *median* over all called SNPs stays at 2-4 even for libraries that cover highly expressed
genes very deeply. The distribution is strongly right-skewed and its median is the wrong summary of
it. Across the 53 panNET RNA-seq libraries used to calibrate this gate, not one sample reached a
median DP of 5, so any non-trivial floor fails every sample. The value is still computed and written
to `qc/sample_qc.tsv` for inspection. **Raise it (e.g. `10`) for WES/WGS input**, where median depth
genuinely is a quality measure. `het_range`'s upper bound is `0.70` rather than `0.60` for the same
reason: low-depth heterozygote over-calling plus allele-specific expression push good RNA-seq
libraries to ~0.6 (observed range across those 53 samples: 0.154-0.609).

**Variant-level QC** is unchanged: study call rate (`params.max_missing`), per-superpopulation
reference MAF (`params.ref_maf`), study coverage (`params.study_coverage`), long-range-LD exclusion,
and LD pruning. Pooled MAF/HWE on the merged panel stay off by default because on a multi-ancestry
panel they preferentially discard ancestry-informative markers.

**Strand-ambiguous SNPs.** `params.exclude_ambiguous` (default `true`) drops palindromic A/T and C/G
SNPs from the study/1000G intersection. SHAKIRA is not structurally exposed to strand error - study
and reference genotypes are both called against GRCh38, variants are matched on exact
`CHROM_POS_REF_ALT` IDs assigned before merging, and every PLINK step runs `--keep-allele-order`, so
no strand or A1/A2 flip is ever attempted and a mismatched variant simply fails to intersect - but
excluding these SNPs removes the whole class of risk at negligible cost. The dropped IDs are written
to `merge/overlap_snps.ambiguous_dropped.txt`; set the flag to `false` to retain them.

## Outputs (under `work_root`/`outdir`, default `<work_root>/amr`)

- `qc/sample_qc.tsv` — per-sample input QC metrics with PASS/FAIL and the reason for each exclusion
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
           get_example_data.sh   downloads the public CCLE example dataset
workflow/  Snakefile  scripts/    the workflow + its helper scripts
           scripts/sample_qc_one.sh  sample_qc_gate.py  drop_ambiguous.sh
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
