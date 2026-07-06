# Preflight checklist

Companion to the top-level `README.md`. Confirms what is configured and what you must provide
before the first run. Reference/param single source of truth: `config.yaml`.

## Configured (portable)

- **`config.yaml`** — reference paths (built from `path_anchors`; 1KG panel + index REQUIRED),
  locked analysis params, caller toggle, per-caller VCF glob, conda-env map. No hardcoded
  site paths: set the three anchors (or export `ANCESTRY_WORK` / `ANCESTRY_PANEL` / `ANCESTRY_GATK`).
- **Conda envs** — `environment.yml` (one standalone controller env) plus per-rule `envs/*.yaml`
  (bcftools, plink, admixture, r-tidyverse, python, bedtools, gatk, nextflow).
- **Execution profiles** — `profiles/slurm/` (cluster; partitions/account set here, not in the
  Snakefile) and `profiles/local/` (single machine, no scheduler).
- **`refs/prep_refs.sh`** — (1) hard-requires the 1KG panel + index, (2) renames the GATK
  known-sites from chr-prefixed to bare contigs + reheaders to the FASTA, (3) builds/validates
  FASTA `.fai`/`.dict` and the panel table. Idempotent. Reads the same env anchors.
- **`TOOLS_VERSIONS.txt`** — as-run vs standardized conda pins; `refs/capture_versions.sh` records
  exact installed builds.
- **Upstream calling** — `Validation/run_rnavar.slurm` / `run_sarek.slurm` (nf-core rnavar 1.3.0 /
  sarek 3.9.0) + `nfcore_hpc.config` (dynamic memory, OOM-retry, resource limits; generalized to
  both callers; partition via `NFCORE_QUEUE`). `map_srr_to_ancestry.py` builds the matching sheet.
- **`init_nfcore_singularity.sh`** — one-time bootstrap: nextflow → nf-core → pull rnavar & sarek →
  pre-stage Singularity images.

## Provide before the first run

1. **Set the anchors** in `config.yaml` `path_anchors` (or export `ANCESTRY_WORK`,
   `ANCESTRY_PANEL`, `ANCESTRY_GATK`). See the README "Configure" section.
2. **1000G panel** `1kgenome.vcf.gz` + its `.tbi`/`.csi` under `panel_dir`. Required-panel gate in
   `prep_refs.sh`. (Provenance in `config.yaml`: 1000G 20190312 biallelic release, per-chr concat.)
3. **Ensembl FASTA** (`Homo_sapiens.GRCh38.dna.primary_assembly.fa`) and **GTF v112** under
   `work_root/ref/...`; run `refs/prep_refs.sh` to build `.fai`/`.dict` and the GATK `_rename` sites.
4. **One example caller VCF** to confirm the glob: rnavar `*.haplotypecaller.filtered.vcf.gz` vs
   sarek `*.haplotypecaller.vcf.gz` (confirm sarek's exact 3.9.0 output name for your run).
5. **Bootstrap nf-core** with `init_nfcore_singularity.sh` before submitting the caller job.

Dry-run to preview the DAG once anchors are set: `snakemake -s workflow/Snakefile -n`.

## Open decisions

- **Caller default**: `rnavar` (STAR + SplitNCigarReads; correct for CCLE RNA-seq) vs `sarek` (DNA
  germline; off-label for RNA). Kept `rnavar`.
- **K3 high-quality-reference toggle**: the optional `<98%` exclusion list is exposed as
  `reference.kg_exclude` (default OFF; admixed samples already dropped by `exclude_subpops` +
  per-superpop MAF). Legacy builder lives under `K3/` (provenance).

## Repo layout

```
README.md  config.yaml  environment.yml  TOOLS_VERSIONS.txt  PREFLIGHT_snakemake.md
run_snakemake.sh  init_nfcore_singularity.sh
envs/       bcftools|plink|admixture|r-tidyverse|python|bedtools|gatk|nextflow .yaml
profiles/   slurm/config.yaml  local/config.yaml
refs/       prep_refs.sh  capture_versions.sh
resources/  1K_pops.txt  high-LD-regions-hg38-GRCh38.txt        # committed reference tables
workflow/   Snakefile
            scripts/ ingest_one.sh ref_maf_per_superpop.sh pop.py analyze_ancestry.R plot_cv_error.R
Validation/ run_rnavar.slurm  run_sarek.slurm  rnavar_params.yaml  sarek_params.yaml
            nfcore_hpc.config  map_srr_to_ancestry.py  ccle_ancestry_master.csv  archive/
cheaha/     legacy standalone driver (provenance; see cheaha/README.md)  archive/
K3/         legacy high-quality-reference scripts (provenance; see K3/README.md)
```
