"""Shared loading/processing of C model output for results.py (figures) and tables.py (tables).

Command-line options (parsed by both scripts):
  --cdir DIR    directory with the model's vars*.csv output (default ../c/output)
  --out DIR     where figures/tables are written (default output)
  --ghh         use the endogenous-labor (GHH) sensitivity runs, i.e. files with the _ghh suffix;
                output names get the same suffix so they do not overwrite the baseline results
  --pz lo|hi    use the sourcing-friction sensitivity runs solved with -z lo or -z hi (files with
                the _phiZ_lo / _phiZ_hi suffix; PHIZ below); outputs get the same suffix. A number
                instead selects a run solved with -z <number> (suffix _pz<number>). Combines with --ghh.
  --partial     allow expected model runs to be missing. Only into a fresh --out directory: the
                skipped inputs and outputs are listed in PARTIAL_MANIFEST.txt there, so a partial
                set can never sit next to stale figures/tables from an earlier complete run.

Without --partial, every expected run must be present (all scenarios for the baseline; the
headline set SENS_RUNS with --ghh or --pz), otherwise the scripts stop before writing anything.
Model output must come from the revised C code: files lacking REQUIRED_COLS (written before the
GDP/consumption/NFA columns were added) are rejected rather than reinterpreted.
"""
import argparse
import os
import numpy as np
import pandas as pd

SLABELS = ['Oil', 'Steel', 'Toys', 'Cars']
SCENARIOS = ['Oil', 'Steel', 'Toys', 'Cars', 'All']
SCEN_CODES = [0, 1, 3, 4, 6]          # -s codes for the five single/all-goods policies

# Columns written by the revised model code (eba2c21 and later). Older files define ngdp on the
# production side and have no CES consumption aggregate, so they are not comparable.
REQUIRED_COLS = ['ngdp_prod', 'b', 'cagg', 'pc_exact', 'welfare']

# Runs that make up each sensitivity set (GHH labor, alternative phi_Z): the five single/all-goods
# policies and the 2025 schedule. Keys as in load_results/_trump2025.
SENS_RUNS = {(25, s, 2, 0, 0, 4) for s in SCEN_CODES} | {('trump2025', 0, 0, 4)}

MANIFEST = 'PARTIAL_MANIFEST.txt'
_missing = []        # expected model files not found, recorded by the loaders


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cdir', default='../c/output')
    ap.add_argument('--out', default='output')
    ap.add_argument('--ghh', action='store_true')
    ap.add_argument('--pz', default=None, help='lo, hi, or a phi_Z value')
    ap.add_argument('--partial', action='store_true')
    args = ap.parse_args()
    args.suffix = suffix(args.ghh, args.pz)
    args.sens = args.ghh or args.pz is not None
    os.makedirs(args.out, exist_ok=True)
    if args.partial:
        manifest = os.path.join(args.out, MANIFEST)
        if not os.path.exists(manifest):
            if os.listdir(args.out):
                raise SystemExit('--partial writes only to a fresh (empty) --out directory; %s is not empty'
                                 % args.out)
            with open(manifest, 'w') as f:
                f.write('Partial results: generated from incomplete model output in %s\n' % args.cdir)
                f.write('Items below were skipped; everything else in this directory was generated.\n')
    return args


# Alternative sourcing adjustment costs (calibrate_phiz.py), as PHIZ_LO / PHIZ_HI in the C code's
# globals.h: short-run Armington elasticity 1.0 (lo, weak frictions) and 0.2 (hi, strong frictions).
PHIZ = {'lo': 0.6585, 'hi': 12.44}
PHIZ_BASE = 2.9


def phiz_value(pz=None):
    """phi_Z for a --pz argument: None (baseline), 'lo', 'hi', or a number."""
    return PHIZ_BASE if pz is None else PHIZ[pz] if pz in PHIZ else float(pz)


def suffix(ghh=False, pz=None):
    """File-name suffix used by the C code: _ghh for GHH labor; _phiZ_lo / _phiZ_hi for -z lo / hi,
    _pz<phi> (as "%g") for -z <number>."""
    sz = '' if pz is None else '_phiZ_' + pz if pz in PHIZ else '_pz%g' % float(pz)
    return ('_ghh' if ghh else '') + sz


def _record(args, line):
    if args.partial:
        with open(os.path.join(args.out, MANIFEST), 'a') as f:
            f.write(line + '\n')


def check_inputs(args, tool):
    """Stop unless every expected run was found (or --partial, which records what is missing)."""
    if not _missing:
        return
    if not args.partial:
        raise SystemExit('%s: expected model output missing (nothing written):\n  %s\n'
                         'Run the missing scenarios, or use --partial with a fresh --out directory.'
                         % (tool, '\n  '.join(_missing)))
    for p in _missing:
        _record(args, '%s: missing input %s' % (tool, p))


def note_skipped(args, tool, what):
    print('%s: skipping %s' % (tool, what))
    _record(args, '%s: skipped %s' % (tool, what))


def process_results(models):
    """Add derived columns. Raw series are left untouched; `<col>_y` columns hold nx series
    normalized by the same period's nominal GDP."""
    base = models[0]
    for m in models:
        m['nx'] = m['nx1'] + m['nx2']
        m['rer'] = m['rer1'] * base['ex1'].values[0] + m['rer2'] * base['ex2'].values[0]
        m['lg'] = m['l0'] + m['l1'] + m['l2'] + m['l3']
        m['lgc'] = m['lg'] + m['l5']          # goods + construction
        for s in range(5):
            for v in ['nxs', 'exs', 'ims', 'rexs', 'rims']:
                m['%s%d' % (v, s)] = m['%s%d-1' % (v, s)] + m['%s%d-2' % (v, s)]
        m['img'] = m['ims0'] + m['ims1'] + m['ims2'] + m['ims3']
        m['exg'] = m['exs0'] + m['exs1'] + m['exs2'] + m['exs3']
        m['nxg'] = m['nxs0'] + m['nxs1'] + m['nxs2'] + m['nxs3']
        m['nxsvc'] = m['nxs4']
        m['rex'] = m['rex1'] + m['rex2']
        m['rim'] = m['rim1'] + m['rim2']
        m['rexg'] = m['rexs0'] + m['rexs1'] + m['rexs2'] + m['rexs3']
        m['rimg'] = m['rims0'] + m['rims1'] + m['rims2'] + m['rims3']
        m['nc'] = m['pc_exact'] * m['cagg']                # nominal consumption = sum_s p_s c_s
        m['ni'] = m['pi'] * m['ii']                        # nominal investment
        for col in [c for c in m.columns if 'nx' in c and not c.endswith('_y')]:
            m[col + '_y'] = m[col] / m['ngdp']             # ngdp: expenditure-side nominal GDP
        m['b_y'] = m['b'] / m['ngdp']
    return models


def _read(path):
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    lacking = [c for c in REQUIRED_COLS if c not in df.columns]
    if lacking:
        raise SystemExit('%s was written by model code that predates the revision (no %s); '
                         'regenerate it with the current code.' % (path, ', '.join(lacking)))
    return df


def _load_pair(key, base, cf, args):
    b, x = _read(base), _read(cf)
    if b is None or x is None:
        if not args.sens or key in SENS_RUNS:        # expected for this set of results
            _missing.extend(p for p, v in [(base, b), (cf, x)] if v is None and p not in _missing)
        return None
    return process_results([b, x])


def load_results(t, s, c, r, d, a, args):
    """Baseline and counterfactual US paths for a scalar-tariff scenario; None if not run."""
    return _load_pair((t, s, c, r, d, a),
                      os.path.join(args.cdir, 'vars0_usa_a%d%s.csv' % (a, args.suffix)),
                      os.path.join(args.cdir, 'vars1_usa_t%d_s%d_c%d_r%d_d%d_a%d%s.csv'
                                   % (t, s, c, r, d, a, args.suffix)), args)


def load_results_trump2025(r, d, a, args):
    return _load_pair(('trump2025', r, d, a),
                      os.path.join(args.cdir, 'vars0_usa_a%d%s.csv' % (a, args.suffix)),
                      os.path.join(args.cdir, 'vars1_usa_trump2025_r%d_d%d_a%d%s.csv'
                                   % (r, d, a, args.suffix)), args)


def have(*dfs):
    """True if every argument is present. An argument is a loaded scenario ([baseline, counterfactual]
    DataFrames), None (not run), or a list of scenarios."""
    for d in dfs:
        if d is None:
            return False
        if isinstance(d, list) and not isinstance(d[0], pd.DataFrame) and any(x is None for x in d):
            return False
    return True


def pct_chg(df, variable):
    return (100 * (df[1][variable] / df[0][variable] - 1)).values


def diff_pct_lg(df, variable):
    return (100 * (df[1][variable] - df[0][variable]) / df[0]['lg']).values


def pp_chg(df, variable):
    return (100 * (df[1][variable] - df[0][variable])).values
