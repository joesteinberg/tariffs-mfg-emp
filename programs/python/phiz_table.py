"""Sourcing-friction sensitivity table: goods-employment transition under alternative phi_Z.

usage: python3 phiz_table.py [--pz lo --pz hi] [--cdir ../c/output] [--out output]

Rows: the all-goods and 2025 experiments, each at the baseline phi_Z (2.9, files without suffix)
and at every --pz value (default lo and hi: files with suffix _phiZ_lo / _phiZ_hi, solved with
model.bin -z lo / -z hi; the values, from calibrate_phiz.py, are PHIZ in model_io.py). Columns:
  phi_Z            sourcing adjustment cost
  calib_elast      the calibration statistic at this phi_Z: US impact Armington elasticity (`ae`
                   at t = 0) in the symmetric all-goods experiment (vars1_usa_t25_s6_c2_r2_d0_a4*)
  exp_elast        the same statistic in the experiment of this row (CSV only: with tariffs that
                   differ across partners or sectors, as in 2025, some relative prices barely move
                   and the trade-weighted ratio is erratic)
  impact, trough, trough_yr, lr
                   goods employment, % change from the no-tariff baseline (lr = final period)
  recovery_yr      first year from which goods employment stays above baseline ('--' if never)
  conv_yr          years to converge to the new long run: first t from which
                   |g_t - g_LR| <= 10% of max_t |g_t - g_LR| for all later t
  last_step        |g_NT - g_NT-1|, to confirm the path has converged within the horizon (CSV only)
Writes tab_phiz.csv and tab_phiz.tex.
"""
import argparse
import os
import numpy as np
import pandas as pd
from model_io import _load_pair, suffix, phiz_value, PHIZ_BASE, recovery_year, conv_year

ap = argparse.ArgumentParser()
ap.add_argument('--cdir', default='../c/output')
ap.add_argument('--out', default='output')
ap.add_argument('--pz', action='append', help='lo, hi, or a phi_Z value (default: lo and hi)')
args = ap.parse_args()
args.pz = args.pz or ['lo', 'hi']
os.makedirs(args.out, exist_ok=True)


class _A:                       # minimal args object for model_io's loader
    ghh = False
    sens = True


def load(exp, pz):
    a = _A()
    sfx = suffix(False, pz)
    base = os.path.join(args.cdir, 'vars0_usa_a4%s.csv' % sfx)
    if exp == 'All goods':
        cf = 'vars1_usa_t25_s6_c2_r0_d0_a4%s.csv' % sfx
        key = (25, 6, 2, 0, 0, 4)
    else:
        cf = 'vars1_usa_trump2025_r0_d0_a4%s.csv' % sfx
        key = ('trump2025', 0, 0, 4)
    df = _load_pair(key, base, os.path.join(args.cdir, cf), a)
    if df is None:
        raise SystemExit('phiz_table.py: missing %s or %s' % (base, cf))
    return df


rows = []
for exp in ['All goods', '2025']:
    for pz in sorted([None] + args.pz, key=phiz_value):   # weak to strong
        b, x = load(exp, pz)
        g = 100 * (x['lg'] / b['lg'] - 1).values
        trough = g.min()
        sym = os.path.join(args.cdir, 'vars1_usa_t25_s6_c2_r2_d0_a4%s.csv' % suffix(False, pz))
        if not os.path.exists(sym):
            raise SystemExit('phiz_table.py: missing calibration run %s' % sym)
        rows.append({'experiment': exp, 'phi_Z': phiz_value(pz),
                     'calib_elast': pd.read_csv(sym)['ae'].values[0],
                     'exp_elast': x['ae'].values[0],
                     'impact': g[0],
                     'trough': trough if trough < 0 else np.nan,
                     'trough_yr': int(g.argmin()) if trough < 0 else np.nan,
                     'recovery_yr': recovery_year(g),
                     'conv_yr': conv_year(g),
                     'lr': g[-1],
                     'last_step': abs(g[-1] - g[-2])})
tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(args.out, 'tab_phiz.csv'), index=False)

def frictions(r):
    """Row label: strength of the frictions and the short-run elasticity phi_Z was set to hit."""
    lab = 'Baseline' if r['phi_Z'] == PHIZ_BASE else 'Weak' if r['phi_Z'] < PHIZ_BASE else 'Strong'
    return '%s (%.2f)' % (lab, r['calib_elast'])


head = ['Experiment', 'Frictions (SR elast.)', '$\\phi_Z$', 'Impact', 'Trough', 'Trough yr',
        'Recovery yr', 'Conv. yr', 'Long run']
lines = ['\\begin{tabular}{llrrrrrrr}', '\\toprule', ' & '.join(head) + ' \\\\', '\\midrule']
for i, r in tab.iterrows():
    if i and r['experiment'] != tab['experiment'][i - 1]:
        lines.append('\\midrule')
    cell = lambda v, f: '--' if isinstance(v, float) and np.isnan(v) else f % v
    lines.append(' & '.join([r['experiment'], frictions(r), '%g' % r['phi_Z'],
                             '%.2f' % r['impact'], cell(r['trough'], '%.2f'), cell(r['trough_yr'], '%d'),
                             cell(r['recovery_yr'], '%d'), '%d' % r['conv_yr'], '%.2f' % r['lr']]) + ' \\\\')
lines += ['\\bottomrule', '\\end{tabular}']
with open(os.path.join(args.out, 'tab_phiz.tex'), 'w') as f:
    f.write('\n'.join(lines) + '\n')
print(tab.round(3).to_string(index=False))
