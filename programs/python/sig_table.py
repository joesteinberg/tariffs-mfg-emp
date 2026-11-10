"""Long-run elasticity sensitivity table: results under alternative Armington elasticity vectors.

usage: python3 sig_table.py [--sig NAME ...] [--cdir ../c/output] [--out output]

Rows: the baseline vector (files without suffix) and each --sig vector (default: cp, krs, half, u2,
umean, u5; files with suffix _sig_NAME, solved with model.bin -e NAME, which uses the phi_Z
re-targeted for that vector; see the C code's globals.h). Each vector's phi_Z is set so that the
short-run elasticity in the symmetric calibration experiment equals the baseline's.

Columns (all relative to the no-tariff equilibrium of the same vector; LR = final period):
  sig_oil..sig_cars   long-run elasticities for intermediate / final sourcing
  phiz, sr_elast      sourcing adjustment cost and achieved short-run elasticity (`ae` at t = 0 in
                      vars1_usa_t25_s6_c2_r2_d0_a4*), the calibration target
  lr_elast            trade-weighted long-run elasticity (`ae` at t = NT in the same run)
  emp_<policy>        LR goods employment, % change
  gdp_<policy>        LR real GDP, % change
  imp_, trough_, trough_yr_, rec_yr_, conv_yr_ <policy>
                      goods-employment dynamics, defined as in phiz_table.py: impact, minimum and
                      its year (if negative), first year from which it stays above baseline, and
                      years to converge to within 10% of its maximum deviation from the long run
  exp_elast_<policy>  short-run elasticity (`ae` at t = 0) in that experiment (CSV only: with tariffs
                      that differ across partners or sectors, some relative prices barely move and
                      the trade-weighted ratio is erratic)
  ranking             single-sector tariffs ordered by LR goods employment, best first
Writes tab_sig.csv (everything), tab_sig.tex (LR goods employment and real GDP), and
tab_sig_dyn.tex (goods-employment transition under tariffs on all goods and the 2025 schedule).
"""
import argparse
import os
import numpy as np
import pandas as pd
from model_io import _load_pair, suffix, SCEN_CODES, recovery_year, conv_year

VECTORS = {   # name: (label, intermediate, final), as in globals.h
    'base': ('Baseline', [8, 2, 4, 1.05], [8, 2, 4, 1.05]),
    'cp': ('Caliendo--Parro', [20.13, 3.87, 12.82, 3.14], [20.13, 3.87, 12.82, 3.14]),
    'krs': ('KRS (3 int., 2 final)', [3] * 4, [2] * 4),
    'half': ('Half dispersion', [5.50, 2.75, 3.89, 1.99], [5.50, 2.75, 3.89, 1.99]),
    'u2': ('Uniform 2', [2] * 4, [2] * 4),
    'umean': ('Uniform 3.25', [3.25] * 4, [3.25] * 4),
    'u5': ('Uniform 5', [5] * 4, [5] * 4),
}
SINGLE = ['Oil', 'Steel', 'Toys', 'Cars']
POLICIES = SINGLE + ['All', '2025']

ap = argparse.ArgumentParser()
ap.add_argument('--cdir', default='../c/output')
ap.add_argument('--out', default='output')
ap.add_argument('--sig', action='append', help='vector name (default: all alternatives)')
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)
names = ['base'] + (args.sig or [v for v in VECTORS if v != 'base'])


class _A:                       # minimal args object for model_io's loader
    ghh = False
    sens = True


def load(policy, sig):
    sfx = suffix(False, None, sig)
    base = os.path.join(args.cdir, 'vars0_usa_a4%s.csv' % sfx)
    if policy == '2025':
        cf, key = 'vars1_usa_trump2025_r0_d0_a4%s.csv' % sfx, ('trump2025', 0, 0, 4)
    else:
        s = dict(zip(['Oil', 'Steel', 'Toys', 'Cars', 'All'], SCEN_CODES))[policy]
        cf, key = 'vars1_usa_t25_s%d_c2_r0_d0_a4%s.csv' % (s, sfx), (25, s, 2, 0, 0, 4)
    df = _load_pair(key, base, os.path.join(args.cdir, cf), _A())
    if df is None:
        raise SystemExit('sig_table.py: missing %s or %s' % (base, cf))
    return df


def pct(b, x, v):
    return 100 * (x[v] / b[v] - 1).values


def phiz_of(sig):
    """phi_Z the run used, from its log header (the C code prints it)."""
    if sig == 'base':
        return 2.9
    log = os.path.join(args.cdir, 'log_t25_s6_c2_r0_d0_a4%s.txt' % suffix(False, None, sig))
    with open(log) as f:
        for line in f:
            if 'Sourcing adjustment cost phi_Z' in line:
                return float(line.split('=')[1].split()[0])
    raise SystemExit('sig_table.py: no phi_Z line in %s' % log)


rows = []
for sig in names:
    label, si, sf = VECTORS[sig]
    sym = pd.read_csv(os.path.join(args.cdir, 'vars1_usa_t25_s6_c2_r2_d0_a4%s.csv' % suffix(False, None, sig)))
    r = {'vector': sig, 'label': label, 'phiz': phiz_of(sig),
         'sr_elast': sym['ae'].values[0], 'lr_elast': sym['ae'].values[-1]}
    for k, nm in enumerate(['oil', 'steel', 'toys', 'cars']):
        r['sig_%s' % nm] = '%g' % si[k] if si[k] == sf[k] else '%g/%g' % (si[k], sf[k])
    for p in POLICIES:
        b, x = load(p, sig)
        g = pct(b, x, 'lg')
        r['emp_%s' % p] = g[-1]
        r['gdp_%s' % p] = pct(b, x, 'rgdp')[-1]
        r['exp_elast_%s' % p] = x['ae'].values[0]
        r['imp_%s' % p] = g[0]
        r['trough_%s' % p] = g.min() if g.min() < 0 else np.nan
        r['trough_yr_%s' % p] = int(g.argmin()) if g.min() < 0 else np.nan
        r['rec_yr_%s' % p] = recovery_year(g)
        r['conv_yr_%s' % p] = conv_year(g)
    r['ranking'] = ' > '.join(sorted(SINGLE, key=lambda p: -r['emp_%s' % p]))
    rows.append(r)
tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(args.out, 'tab_sig.csv'), index=False)

ncol = 1 + 2 * len(POLICIES)
lines = ['\\begin{tabular}{l' + 'r' * (ncol - 1) + '}', '\\toprule',
         ' & \\multicolumn{%d}{c}{LR goods employment (\\%%)} & \\multicolumn{%d}{c}{LR real GDP (\\%%)} \\\\'
         % (len(POLICIES), len(POLICIES)),
         '\\cmidrule(lr){2-%d}\\cmidrule(lr){%d-%d}' % (1 + len(POLICIES), 2 + len(POLICIES), ncol),
         'Elasticities & ' + ' & '.join(POLICIES * 2) + ' \\\\', '\\midrule']
for _, r in tab.iterrows():
    lines.append(' & '.join([r['label']] + ['%.2f' % r['emp_%s' % p] for p in POLICIES]
                            + ['%.2f' % r['gdp_%s' % p] for p in POLICIES]) + ' \\\\')
lines += ['\\bottomrule', '\\end{tabular}']
with open(os.path.join(args.out, 'tab_sig.tex'), 'w') as f:
    f.write('\n'.join(lines) + '\n')

cell = lambda v, f: '--' if isinstance(v, float) and np.isnan(v) else f % v
head = ['Experiment', 'Elasticities', '$\\phi_Z$', 'Impact', 'Trough', 'Trough yr',
        'Recovery yr', 'Conv. yr', 'Long run']
lines = ['\\begin{tabular}{llrrrrrrr}', '\\toprule', ' & '.join(head) + ' \\\\', '\\midrule']
for k, (p, lab) in enumerate([('All', 'All goods'), ('2025', '2025')]):
    if k:
        lines.append('\\midrule')
    for _, r in tab.iterrows():
        lines.append(' & '.join([lab, r['label'], '%g' % r['phiz'], '%.2f' % r['imp_%s' % p], cell(r['trough_%s' % p], '%.2f'),
                                 cell(r['trough_yr_%s' % p], '%d'), cell(r['rec_yr_%s' % p], '%d'),
                                 '%d' % r['conv_yr_%s' % p], '%.2f' % r['emp_%s' % p]]) + ' \\\\')
lines += ['\\bottomrule', '\\end{tabular}']
with open(os.path.join(args.out, 'tab_sig_dyn.tex'), 'w') as f:
    f.write('\n'.join(lines) + '\n')

with pd.option_context('display.width', 200):
    print(tab[['vector', 'phiz', 'sr_elast', 'lr_elast'] + ['emp_%s' % p for p in POLICIES] + ['ranking']]
          .round(3).to_string(index=False))
