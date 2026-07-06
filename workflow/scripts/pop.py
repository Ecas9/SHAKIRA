#!/usr/bin/env python3
"""Build an ADMIXTURE .pop label file aligned 1:1 with a PLINK .fam file.

Supervised ADMIXTURE requires a `.pop` file with EXACTLY one line per individual
in the .fam, in the SAME order. Each line holds that individual's (super)population
label, or a blank line for individuals to be estimated (e.g. the study samples).

This replaces the previous version, which:
  * iterated a set() of sample IDs, destroying the .fam row order (every label
    landed on the wrong individual);
  * never skipped the header row of the sample->population table; and
  * hardcoded all input/output paths.

Usage (defaults assume you run it from step_4/):
    python3 pop.py                       # prunedData.fam -> prunedData.pop
    python3 pop.py --fam other.fam --out other.pop
    python3 pop.py --label pop           # write sub-population codes instead of super
"""
import argparse
import os
import sys

# Fallback 1000 Genomes population -> super-population map. Prefer --pop-table
# (e.g. the committed 20131219.populations.tsv) so this list never drifts.
POP_TO_SUPER = {
    'CHB': 'EAS', 'JPT': 'EAS', 'CHS': 'EAS', 'CDX': 'EAS', 'KHV': 'EAS', 'CHD': 'EAS',
    'CEU': 'EUR', 'TSI': 'EUR', 'GBR': 'EUR', 'FIN': 'EUR', 'IBS': 'EUR',
    'YRI': 'AFR', 'LWK': 'AFR', 'GWD': 'AFR', 'MSL': 'AFR', 'ESN': 'AFR',
    'ASW': 'AFR', 'ACB': 'AFR',
    'MXL': 'AMR', 'PUR': 'AMR', 'CLM': 'AMR', 'PEL': 'AMR',
    'GIH': 'SAS', 'PJL': 'SAS', 'BEB': 'SAS', 'STU': 'SAS', 'ITU': 'SAS',
}


def load_pop_to_super(path):
    """Build {population_code: super_population} from a populations TSV.

    Expects the 20131219.populations.tsv layout: column 2 = 'Population Code',
    column 3 = 'Super Population'. The header row is skipped.
    """
    mapping = {}
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3:
                continue
            pop, sup = cols[1].strip(), cols[2].strip()
            if pop and sup:
                mapping[pop] = sup
    return mapping


def load_sample_to_pop(path, skip_header=True):
    """Build {sample_id: population_code} from a 2-column sample->pop table.

    Tolerates either tab- or whitespace-delimited input.
    """
    mapping = {}
    with open(path) as fh:
        if skip_header:
            fh.readline()
        for line in fh:
            line = line.rstrip('\n')
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


def read_fam_order(path):
    """Return the list of family IDs (.fam column 1) in file order."""
    ids = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            ids.append(line.split()[0])
    return ids


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fam', default='prunedData.fam',
                    help='PLINK .fam whose individual order the .pop must match')
    ap.add_argument('--sample-pop', default='1K_pops.txt',
                    help='2-column sample_id<TAB>population_code table')
    ap.add_argument('--no-sample-header', action='store_true',
                    help='set if --sample-pop has no header row to skip')
    ap.add_argument('--pop-table', default=None,
                    help='populations TSV to derive population->super-population '
                         '(e.g. 20131219.populations.tsv); falls back to a built-in map')
    ap.add_argument('--label', choices=['super', 'pop'], default='super',
                    help="write super-population codes (default) or raw population codes")
    ap.add_argument('--out', default=None,
                    help='output .pop path (default: <fam basename>.pop)')
    args = ap.parse_args(argv)

    # Start from the built-in map (covers all 28 1KG population codes, incl. CHD) and
    # overlay the authoritative --pop-table on top. This way a population code present
    # in the sample panel but missing from the table (e.g. CHD, which is not in
    # 20131219.populations.tsv) still resolves to its super-population instead of being
    # left blank — while the table remains the source of truth where it has an entry.
    pop_to_super = dict(POP_TO_SUPER)
    if args.pop_table and os.path.exists(args.pop_table):
        pop_to_super.update(load_pop_to_super(args.pop_table))
    elif args.pop_table:
        sys.stderr.write(f"warning: --pop-table {args.pop_table} not found; "
                         "using built-in population map\n")

    sample_to_pop = load_sample_to_pop(args.sample_pop,
                                       skip_header=not args.no_sample_header)
    fam_ids = read_fam_order(args.fam)
    out_path = args.out or (os.path.splitext(args.fam)[0] + '.pop')

    labeled = 0
    unknown_pops = set()
    with open(out_path, 'w', newline='\n') as out:  # force LF for the Linux ADMIXTURE reader
        for sid in fam_ids:
            pop = sample_to_pop.get(sid)
            if not pop:
                label = ''
            elif args.label == 'pop':
                label = pop
            else:
                label = pop_to_super.get(pop, '')
                if not label:
                    unknown_pops.add(pop)
            if label:
                labeled += 1
            out.write(label + '\n')

    sys.stderr.write(
        f"{out_path}: {len(fam_ids)} individuals, {labeled} labeled, "
        f"{len(fam_ids) - labeled} blank\n")
    if unknown_pops:
        sys.stderr.write("warning: population codes with no super-population "
                         f"mapping (left blank): {sorted(unknown_pops)}\n")
    if labeled == 0:
        sys.stderr.write(
            "ERROR: no individuals were labeled. Check that --sample-pop IDs use "
            "the same namespace as the .fam IDs (e.g. NA*/HG* vs numeric).\n")
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
