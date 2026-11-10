"""Portfolio-adjustment-cost sensitivity: long-run results when net foreign assets are pinned down.

usage: python3 psib_table.py [--psib 0.02 --psib 0.01 --psib 0.005] [--nt 200] [--cdir ../c/output]
                             [--out output]

In the baseline model the long-run equilibrium depends on the net foreign assets (NFA) accumulated
along the transition. With a portfolio adjustment cost psi_B (C option -b, eqm.h) every country's
bonds return to their initial level b0, so the long run is the unique steady state with b = b0,
whatever psi_B; psi_B only governs how fast the economy gets there.

Policies: 25% tariffs on Oil, Steel, Toys, Cars, and all goods (China and the rest of the world),
and the 2025 schedule. Models (columns of the LaTeX table):
  Baseline          files without suffix (endogenous NFA); LR = final period, t = 100
  psi_B = <psi>     files with suffix _psib<psi>_nt<nt>, solved with bin/model_nt<nt>.bin -b <psi>;
                    LR = final period, t = nt
  Steady state      ss1_*_psib<psi>_nt<nt>.csv: the tariff steady state with b = b0, solved directly
                    by the C code (identical for every psi_B; taken from the first --psib and checked
                    against the others)
Statistics (all relative to the no-tariff equilibrium of the same model):
  emp_impact, emp_trough, emp_trough_yr, emp_lr   goods employment, % change
  tb_lr       change in the trade balance, p.p. of GDP (nominal NX / nominal GDP)
  nfa_lr      change in bond holdings, % of no-tariff nominal GDP (CSV only: a check that bonds
              return to b0, where it is zero; the level of b0 itself is not a measured NIIP)
  gdp_lr, cons_lr                                  real GDP and real consumption, % change
  l_<sector>  LR change in sector employment, % of baseline goods employment (goods sectors; these
              add up to emp_lr), and % change for services and construction
Writes tab_psib.csv (everything, long format) and tab_psib.tex (panels: LR goods employment and
LR trade balance).
"""
import argparse
import os
import numpy as np
import pandas as pd
from model_io import _read, process_results, suffix, SLABELS, SCEN_CODES

ap = argparse.ArgumentParser()
ap.add_argument('--cdir', default='../c/output')
ap.add_argument('--out', default='output')
ap.add_argument('--psib', action='append', help='psi_B values (default 0.02, 0.01, 0.005)')
ap.add_argument('--nt', type=int, default=200)
args = ap.parse_args()
args.psib = args.psib or ['0.02', '0.01', '0.005']
os.makedirs(args.out, exist_ok=True)

POLICIES = [(n, 'vars1_usa_t25_s%d_c2_r0_d0_a4%%s.csv' % s) for n, s in zip(SLABELS + ['All goods'], SCEN_CODES)]
POLICIES.append(('2025', 'vars1_usa_trump2025_r0_d0_a4%s.csv'))


def need(path):
    df = _read(path)
    if df is None:
        raise SystemExit('psib_table.py: missing %s' % path)
    return df


def stats(b, x):
    """Statistics for a [baseline, counterfactual] pair (already processed)."""
    g = 100 * (x['lg'] / b['lg'] - 1).values
    trough = g.min()
    r = {'emp_impact': g[0], 'emp_trough': trough if trough < 0 else np.nan,
         'emp_trough_yr': int(g.argmin()) if trough < 0 else np.nan, 'emp_lr': g[-1],
         'tb_lr': 100 * (x['nx_y'].values[-1] - b['nx_y'].values[-1]),
         'nfa_lr': 100 * (x['b'].values[-1] - b['b'].values[-1]) / b['ngdp'].values[-1],
         'gdp_lr': 100 * (x['rgdp'].values[-1] / b['rgdp'].values[-1] - 1),
         'cons_lr': 100 * (x['cagg'].values[-1] / b['cagg'].values[-1] - 1)}
    for s, lab in enumerate(SLABELS):
        r['l_' + lab] = 100 * (x['l%d' % s].values[-1] - b['l%d' % s].values[-1]) / b['lg'].values[-1]
    r['l_Svcs'] = 100 * (x['l4'].values[-1] / b['l4'].values[-1] - 1)
    r['l_Const'] = 100 * (x['l5'].values[-1] / b['l5'].values[-1] - 1)
    return r


rows = []
models = [('Baseline', '', None)] + [('psi_B = %s' % p, suffix(psib=p, nt=args.nt), p) for p in args.psib]
for label, sfx, p in models:
    b0 = need(os.path.join(args.cdir, 'vars0_usa_a4%s.csv' % sfx))
    for pol, pat in POLICIES:
        b, x = process_results([b0.copy(), need(os.path.join(args.cdir, pat % sfx))])
        rows.append(dict(model=label, psib=float(p) if p else 0.0, policy=pol, **stats(b, x)))

# steady state with b = b0 (one row, period NT), the same for every psi_B
ss = {}
for p in args.psib:
    sfx = suffix(psib=p, nt=args.nt)
    b0 = need(os.path.join(args.cdir, 'vars0_usa_a4%s.csv' % sfx)).iloc[[-1]].reset_index(drop=True)
    for pol, pat in POLICIES:
        x = need(os.path.join(args.cdir, ('ss1' + pat[5:]) % sfx))
        b, x = process_results([b0.copy(), x])
        ss.setdefault(pol, []).append(stats(b, x))
for pol, _ in POLICIES:
    r = ss[pol][0]
    for other in ss[pol][1:]:
        dev = max(abs(other[k] - r[k]) for k in r if not (isinstance(r[k], float) and np.isnan(r[k])))
        if dev > 1e-6:
            raise SystemExit('psib_table.py: steady states at b0 differ across psi_B for %s (%g)' % (pol, dev))
    r = dict(r, emp_impact=np.nan, emp_trough=np.nan, emp_trough_yr=np.nan)
    rows.append(dict(model='Steady state', psib=np.nan, policy=pol, **r))

tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(args.out, 'tab_psib.csv'), index=False)

# LaTeX: policies in rows; baseline, each psi_B (final period), and the steady state in columns
order = [m[0] for m in models] + ['Steady state']
head = ['', 'Baseline'] + ['$\\psi_B = %s$' % p for p in args.psib] + ['Steady state']
lines = ['\\begin{tabular}{l' + 'r' * (len(order)) + '}', '\\toprule', ' & '.join(head) + ' \\\\']
panels = [('emp_lr', 'Goods employment (\\% change)'),
          ('tb_lr', 'Trade balance (change, p.p.\\ of GDP)')]
for k, (var, title) in enumerate(panels):
    lines += ['\\midrule', '\\multicolumn{%d}{l}{\\textit{%s. %s}} \\\\' % (len(order) + 1, 'ABC'[k], title)]
    for pol, _ in POLICIES:
        cells = []
        for m in order:
            v = tab[(tab.model == m) & (tab.policy == pol)][var].values[0]
            cells.append('%.2f' % (0.0 if abs(v) < 0.005 else v))
        lines.append(pol + ' & ' + ' & '.join(cells) + ' \\\\')
lines += ['\\bottomrule', '\\end{tabular}']
with open(os.path.join(args.out, 'tab_psib.tex'), 'w') as f:
    f.write('\n'.join(lines) + '\n')

show = tab.pivot(index='policy', columns='model', values='emp_lr')[order].loc[[p for p, _ in POLICIES]]
print('LR goods employment (% change)')
print(show.round(2).to_string())
print(tab.round(3).to_string(index=False))
