#!/usr/bin/env python3
"""Build an ADMIXTURE .pop label file for the K3 high-quality 1KG panel run.

Supervised ADMIXTURE reads ``<prefix>.pop`` as EXACTLY one label per individual in
the companion ``.fam``, in the SAME order. Each line holds that individual's reference
(super-)population label, or an "unknown" sentinel ("-" by default) for individuals
whose ancestry ADMIXTURE should ESTIMATE rather than hold fixed.

This is the K3/JHU panel variant of ``cheaha/pop.py``. It deliberately differs from
``cheaha/pop.py`` in exactly two ways, preserved here as defaults:
  * the unknown sentinel is ``-`` (``cheaha/pop.py`` writes a blank line); and
  * it keys the ``.fam`` on the IID (column 2). Under the pipeline's ``--double-id``
    beds FID==IID, so this is equivalent to ``cheaha/pop.py``'s FID keying.
Everything else (pure stdlib, .fam-order iteration, header skip, guards) mirrors
``cheaha/pop.py``.

This is a rewrite of the original hardcoded-path script, which:
  * fed the population code straight into a 26-entry dict that was MISSING 'CHD', so
    it raised ``KeyError`` and crashed on any CHD individual (1K_pops.txt has CHD);
  * called ``pd.read_csv(path, '\\s+', ...)`` with the separator POSITIONAL, a
    ``TypeError`` on pandas >= 2 (``sep`` is keyword-only); and
  * hardcoded every input/output path (and pulled in pandas + tabulate just to do a
    per-sample DataFrame scan, i.e. an O(n*m) lookup).

Reproduce the original run with:
    python3 pop_jhu.py \\
        --fam /path/to/work/admixture_run/K3_var/prunedData_panel.fam \\
        --out /path/to/work/admixture_run/K3_var/prunedData_panel.pop \\
        --sample-pop /path/to/shakira/resources/1K_pops.txt

For a supervised K=3 run anchored on three continental references only, add e.g.
``--restrict AFR,EUR,EAS`` so AMR/SAS individuals are left as "-" (estimated) instead
of fixed — keeping the number of distinct reference labels equal to K.
"""
import argparse
import os
import sys

# 1000 Genomes population -> super-population, all 28 codes INCLUDING CHD (which the
# original map omitted, causing the KeyError crash). Override/extend via --pop-table.
POP_TO_SUPER = {
    'CHB': 'EAS', 'JPT': 'EAS', 'CHS': 'EAS', 'CDX': 'EAS', 'KHV': 'EAS', 'CHD': 'EAS',
    'CEU': 'EUR', 'TSI': 'EUR', 'GBR': 'EUR', 'FIN': 'EUR', 'IBS': 'EUR',
    'YRI': 'AFR', 'LWK': 'AFR', 'GWD': 'AFR', 'MSL': 'AFR', 'ESN': 'AFR',
    'ASW': 'AFR', 'ACB': 'AFR',
    'MXL': 'AMR', 'PUR': 'AMR', 'CLM': 'AMR', 'PEL': 'AMR',
    'GIH': 'SAS', 'PJL': 'SAS', 'BEB': 'SAS', 'STU': 'SAS', 'ITU': 'SAS',
}


def load_pop_to_super(path):
    """{population_code: super_population} from a 20131219.populations.tsv-style file
    (col 2 = 'Population Code', col 3 = 'Super Population'; header skipped)."""
    mapping = {}
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 3 and cols[1].strip() and cols[2].strip():
                mapping[cols[1].strip()] = cols[2].strip()
    return mapping


def load_sample_to_pop(path, skip_header=True):
    """{sample_id: population_code} from the 2-column SAMPLE_NAME<ws>POPULATION panel.
    Tolerates tab- or whitespace-delimited rows and duplicate rows (last wins)."""
    mapping = {}
    with open(path) as fh:
        if skip_header:
            fh.readline()
        for line in fh:
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                continue
            sample, pop = parts[0].strip(), parts[1].strip()
            if sample:
                mapping[sample] = pop
    return mapping


def read_fam_ids(path, col=1):
    """Individual IDs from a PLINK .fam in file order. col=1 -> IID (2nd field);
    under --double-id FID==IID, so this matches cheaha/pop.py's FID (col=0) keying.
    Falls back to FID for any short/malformed row."""
    ids = []
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if fields:
                ids.append(fields[col] if len(fields) > col else fields[0])
    return ids


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fam', default='prunedData_panel.fam',
                    help='PLINK .fam whose individual order the .pop must match')
    ap.add_argument('--sample-pop', default='1K_pops.txt',
                    help='2-column sample<ws>population panel (SAMPLE_NAME/POPULATION header)')
    ap.add_argument('--no-sample-header', action='store_true',
                    help='set if --sample-pop has no header row to skip')
    ap.add_argument('--pop-table', default=None,
                    help='authoritative populations TSV (e.g. 20131219.populations.tsv) '
                         'overlaid on the built-in map')
    ap.add_argument('--label', choices=['super', 'pop'], default='super',
                    help='write super-population codes (default) or raw population codes')
    ap.add_argument('--restrict', default=None,
                    help='comma-separated super-pops to keep as REFERENCES (e.g. '
                         'AFR,EUR,EAS for a K=3 run); individuals outside the set are '
                         'left as the unknown label (estimated). Requires --label super.')
    ap.add_argument('--unknown-label', default='-',
                    help="sentinel for individuals to ESTIMATE (default '-'; pass '' for "
                         "a blank line, as cheaha/pop.py writes)")
    ap.add_argument('--fam-col', type=int, default=1,
                    help='.fam column to key on (0=FID, 1=IID; equal under --double-id)')
    ap.add_argument('--out', default=None,
                    help='output .pop path (default: <fam basename>.pop)')
    args = ap.parse_args(argv)

    if args.restrict and args.label != 'super':
        ap.error('--restrict filters by super-population and requires --label super')

    # Built-in 28-pop map, with the authoritative --pop-table overlaid where present.
    pop_to_super = dict(POP_TO_SUPER)
    if args.pop_table and os.path.exists(args.pop_table):
        pop_to_super.update(load_pop_to_super(args.pop_table))
    elif args.pop_table:
        sys.stderr.write(f"warning: --pop-table {args.pop_table} not found; "
                         "using built-in population map\n")

    restrict = None
    if args.restrict:
        restrict = {s.strip().upper() for s in args.restrict.split(',') if s.strip()}

    sample_to_pop = load_sample_to_pop(args.sample_pop,
                                       skip_header=not args.no_sample_header)
    fam_ids = read_fam_ids(args.fam, col=args.fam_col)
    out_path = args.out or (os.path.splitext(args.fam)[0] + '.pop')
    unk = args.unknown_label

    labeled = 0
    unknown_pops = set()
    with open(out_path, 'w', newline='\n') as out:  # force LF for the Linux ADMIXTURE reader
        for sid in fam_ids:
            pop = sample_to_pop.get(sid)
            if not pop:
                label = unk
            elif args.label == 'pop':
                label = pop
            else:
                label = pop_to_super.get(pop, '')
                if not label:
                    unknown_pops.add(pop)
                    label = unk
            if restrict is not None and label not in restrict:
                label = unk
            if label and label != unk:
                labeled += 1
            out.write(label + '\n')

    n = len(fam_ids)
    sys.stderr.write(f"{out_path}: {n} individuals, {labeled} labeled, "
                     f"{n - labeled} '{unk or 'blank'}'\n")
    if unknown_pops:
        sys.stderr.write("warning: population codes with no super-population mapping "
                         f"(left unknown): {sorted(unknown_pops)}\n")
    if restrict is not None:
        sys.stderr.write(f"restricted references to {sorted(restrict)} (others estimated)\n")
    if labeled == 0:
        sys.stderr.write("ERROR: no individuals were labeled. Check that --sample-pop IDs "
                         "share the .fam ID namespace (e.g. NA*/HG*).\n")
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
