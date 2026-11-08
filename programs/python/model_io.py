"""Shared loading/processing of C model output for results.py (figures) and tables.py (tables).

Command-line options (parsed by both scripts):
  --cdir DIR    directory with the model's vars*.csv output (default ../c/output)
  --out DIR     where figures/tables are written (default output)
  --ghh         use the endogenous-labor (GHH) sensitivity runs, i.e. files with the _ghh suffix;
                output names get the same suffix so they do not overwrite the baseline results
"""
import argparse
import os
import numpy as np
import pandas as pd

SLABELS = ['Oil', 'Steel', 'Toys', 'Cars']
SCENARIOS = ['Oil', 'Steel', 'Toys', 'Cars', 'All']
SCEN_CODES = [0, 1, 3, 4, 6]          # -s codes for the five single/all-goods policies


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cdir', default='../c/output')
    ap.add_argument('--out', default='output')
    ap.add_argument('--ghh', action='store_true')
    args = ap.parse_args()
    args.suffix = '_ghh' if args.ghh else ''
    os.makedirs(args.out, exist_ok=True)
    return args


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
        if 'cagg' not in m:                  # pre-revision output: fall back to summed goods consumption
            m['cagg'] = m['c']
        if 'pc_exact' in m:
            m['nc'] = m['pc_exact'] * m['cagg']            # nominal consumption = sum_s p_s c_s
        m['ni'] = m['pi'] * m['ii']                        # nominal investment
        for col in [c for c in m.columns if 'nx' in c and not c.endswith('_y')]:
            m[col + '_y'] = m[col] / m['ngdp']
        if 'b' in m:
            m['b_y'] = m['b'] / m['ngdp']
    return models


def _read(path):
    return pd.read_csv(path) if os.path.exists(path) else None


def load_results(t, s, c, r, d, a, args):
    """Baseline and counterfactual US paths for a scalar-tariff scenario; None if not run."""
    b = _read(os.path.join(args.cdir, 'vars0_usa_a%d%s.csv' % (a, args.suffix)))
    x = _read(os.path.join(args.cdir, 'vars1_usa_t%d_s%d_c%d_r%d_d%d_a%d%s.csv' % (t, s, c, r, d, a, args.suffix)))
    if b is None or x is None:
        return None
    return process_results([b, x])


def load_results_trump2025(r, d, a, args):
    b = _read(os.path.join(args.cdir, 'vars0_usa_a%d%s.csv' % (a, args.suffix)))
    x = _read(os.path.join(args.cdir, 'vars1_usa_trump2025_r%d_d%d_a%d%s.csv' % (r, d, a, args.suffix)))
    if b is None or x is None:
        return None
    return process_results([b, x])


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
