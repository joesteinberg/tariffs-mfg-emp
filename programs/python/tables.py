"""Tables and decompositions from the model output (CSV + LaTeX tabular snippets).

usage: python3 tables.py [--cdir ../c/output] [--out output] [--ghh]

All numbers are computed from the vars*.csv files written by the C model. Changes are relative to
the no-tariff baseline in the same period; "LR" is the final period of the transition (t = NT).
Percent changes are simple percent changes (as in results.py's figures) except in the labor
decomposition, which uses log changes so the pieces add up exactly.

Outputs (each with the --ghh suffix when relevant):
  tab_headline.{csv,tex}        goods employment dynamics, aggregates, trade balance (six policies)
  tab_sectors_lr.{csv,tex}      LR goods employment by sector: own % change and contribution
  tab_goods_cns.{csv,tex}       goods + construction employment
  tab_trade_sectors_lr.{csv,tex} LR change in real exports/imports by sector, % of baseline GDP
  tab_va_sectors_lr.{csv,tex}   LR sectoral real value added, labor, capital, VA per worker
  paths_labor_decomp.csv        dlog goods emp = dlog total emp + dlog goods share (paths)
  paths_external.csv            trade balance (total/goods/services/bilateral), real ex/im,
                                NFA, C/Y, I/Y paths
"""
import os
import numpy as np
import pandas as pd
from model_io import (parse_args, load_results, load_results_trump2025, SLABELS, SCEN_CODES)

args = parse_args()

POLICIES = {'Oil': None, 'Steel': None, 'Toys': None, 'Cars': None, 'All': None, '2025': None}
for name, s in zip(['Oil', 'Steel', 'Toys', 'Cars', 'All'], SCEN_CODES):
    POLICIES[name] = load_results(25, s, 2, 0, 0, 4, args)
POLICIES['2025'] = load_results_trump2025(0, 0, 4, args)
missing = [k for k, v in POLICIES.items() if v is None]
if missing:
    print('tables.py: no output for %s in %s; those rows are omitted' % (', '.join(missing), args.cdir))
POLICIES = {k: v for k, v in POLICIES.items() if v is not None}
if not POLICIES:
    raise SystemExit('tables.py: no model output found in %s' % args.cdir)

SECTORS5 = SLABELS + ['Svcs']


def pct(df, v):
    return 100 * (df[1][v] / df[0][v] - 1).values


def tb(m):
    return 100 * m['nx_y'].values


def first_sustained_positive(g):
    """First period from which the series stays strictly positive; NaN if it ends <= 0."""
    if g[-1] <= 0:
        return np.nan
    neg = np.where(g <= 0)[0]
    return 0 if len(neg) == 0 else int(neg.max() + 1)


def path_stats(g):
    trough = g.min()
    return dict(impact=g[0],
                trough=trough if trough < 0 else np.nan,
                trough_yr=int(g.argmin()) if trough < 0 else np.nan,
                recovery_yr=first_sustained_positive(g),
                lr=g[-1])


def out(name):
    return os.path.join(args.out, name + args.suffix)


def write(df, name, fmt, header=None):
    """CSV with full precision; LaTeX tabular (booktabs) with formatted cells."""
    df.to_csv(out(name) + '.csv')
    cols = list(df.columns)
    lines = ['\\begin{tabular}{l' + 'r' * len(cols) + '}', '\\toprule']
    if header:
        lines += header
    else:
        lines.append(' & ' + ' & '.join(cols) + ' \\\\')
    lines.append('\\midrule')
    for idx, row in df.iterrows():
        cells = []
        for c in cols:
            v = row[c]
            f = fmt.get(c, '%.2f') if isinstance(fmt, dict) else fmt
            if isinstance(v, (float, np.floating)) and np.isnan(v):
                cells.append('--')
            elif isinstance(v, str):
                cells.append(v)
            else:
                cells.append(f % v)
        lines.append(str(idx) + ' & ' + ' & '.join(cells) + ' \\\\')
    lines += ['\\bottomrule', '\\end{tabular}']
    with open(out(name) + '.tex', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------------------------
# Headline numbers

rows = {}
for name, m in POLICIES.items():
    b, x = m
    g = pct(m, 'lg')
    r = {('Goods employment', k): v for k, v in path_stats(g).items()}
    r[('Aggregates (LR)', 'total labor')] = pct(m, 'll')[-1]
    r[('Aggregates (LR)', 'real GDP')] = pct(m, 'rgdp')[-1]
    r[('Aggregates (LR)', 'real consumption')] = pct(m, 'cagg')[-1]
    tb0, tb1 = tb(b), tb(x)
    d = tb1 - tb0
    wider = np.where(d < 0)[0]
    r[('Trade balance/GDP', 'baseline')] = tb0[0]
    r[('Trade balance/GDP', 'impact')] = tb1[0]
    r[('Trade balance/GDP', 'max narrowing (pp)')] = d.max()
    r[('Trade balance/GDP', 'yr deficit exceeds baseline')] = int(wider.min()) if len(wider) else np.nan
    r[('Trade balance/GDP', 'LR')] = tb1[-1]
    rows[name] = r
head = pd.DataFrame(rows).T
head.columns = pd.MultiIndex.from_tuples(head.columns)
flat = head.copy()
flat.columns = ['%s: %s' % c for c in head.columns]
fmt = {c: ('%d' if 'yr' in c else '%.2f') for c in flat.columns}   # year columns as integers
hdr = [' & \\multicolumn{5}{c}{Goods employment (\\%)} & \\multicolumn{3}{c}{LR (\\%)} & \\multicolumn{5}{c}{Trade balance / GDP (\\%)} \\\\',
       '\\cmidrule(lr){2-6}\\cmidrule(lr){7-9}\\cmidrule(lr){10-14}',
       ' & Impact & Trough & Year & Recovers & LR & Labor & GDP & Cons. & Base & Impact & Max narrow. & Yr wider & LR \\\\']
write(flat, 'tab_headline', fmt, hdr)

# ---------------------------------------------------------------------------------------------
# LR goods employment by sector: own % change and contribution (pp of baseline goods emp)

rows = {}
for name, m in POLICIES.items():
    b, x = m
    r = {}
    for k, s in enumerate(SLABELS):
        r['%s own %%' % s] = pct(m, 'l%d' % k)[-1]
    for k, s in enumerate(SLABELS):
        r['%s contrib' % s] = 100 * (x['l%d' % k].values[-1] - b['l%d' % k].values[-1]) / b['lg'].values[-1]
    r['Total'] = pct(m, 'lg')[-1]
    r['Gross reallocation'] = sum(abs(r['%s contrib' % s]) for s in SLABELS)
    rows[name] = r
sec = pd.DataFrame(rows).T
hdr = [' & \\multicolumn{4}{c}{Own change (\\%)} & \\multicolumn{4}{c}{Contribution (\\% of goods emp.)} & & \\\\',
       '\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}',
       ' & ' + ' & '.join(SLABELS) + ' & ' + ' & '.join(SLABELS) + ' & Net & Gross \\\\']
write(sec, 'tab_sectors_lr', '%.2f', hdr)

# ---------------------------------------------------------------------------------------------
# Goods + construction

rows = {}
for name, m in POLICIES.items():
    gs = path_stats(pct(m, 'lg'))
    gc = path_stats(pct(m, 'lgc'))
    rows[name] = {'Goods impact': gs['impact'], 'Goods trough': gs['trough'], 'Goods LR': gs['lr'],
                  'Goods+cons impact': gc['impact'], 'Goods+cons trough': gc['trough'], 'Goods+cons LR': gc['lr'],
                  'Cons LR': pct(m, 'l5')[-1]}
gcn = pd.DataFrame(rows).T
hdr = [' & \\multicolumn{3}{c}{Goods} & \\multicolumn{3}{c}{Goods + construction} & Construction \\\\',
       '\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}',
       ' & Impact & Trough & LR & Impact & Trough & LR & LR \\\\']
write(gcn, 'tab_goods_cns', '%.2f', hdr)

# ---------------------------------------------------------------------------------------------
# LR sectoral trade accounting: change in real exports and imports, % of baseline GDP
# (quantity change valued at baseline producer prices of the exporting country)

rows = {}
for name, m in POLICIES.items():
    b, x = m
    y0 = b['ngdp'].values[-1]
    r = {}
    for k, s in enumerate(SECTORS5):
        # US exports of sector k: valued at US baseline producer price py_k
        r['%s ex' % s] = 100 * (x['rexs%d' % k].values[-1] - b['rexs%d' % k].values[-1]) * b['py%d' % k].values[-1] / y0
        # US imports of sector k: valued at baseline import prices (nominal/real at baseline)
        pim = b['ims%d' % k].values[-1] / b['rims%d' % k].values[-1]
        r['%s im' % s] = 100 * (x['rims%d' % k].values[-1] - b['rims%d' % k].values[-1]) * pim / y0
    rows[name] = r
trd = pd.DataFrame(rows).T
hdr = [' & ' + ' & '.join('\\multicolumn{2}{c}{%s}' % s for s in SECTORS5) + ' \\\\',
       ''.join('\\cmidrule(lr){%d-%d}' % (2 + 2 * k, 3 + 2 * k) for k in range(5)),
       ' & ' + ' & '.join(['Ex & Im'] * 5) + ' \\\\']
trd = trd[[c for s in SECTORS5 for c in ('%s ex' % s, '%s im' % s)]]
write(trd, 'tab_trade_sectors_lr', '%.2f', hdr)

# ---------------------------------------------------------------------------------------------
# GDP decomposition inputs: LR % change in sectoral real VA, labor, capital, VA per worker,
# and each sector's contribution to the change in real GDP (pp). The residual is the part of the
# real-GDP change not accounted for by sectoral VA (sourcing adjustment costs and Armington
# aggregation of intermediates).

ALL6 = SLABELS + ['Svcs', 'Cons']
rows = {}
for name, m in POLICIES.items():
    b, x = m
    rg0 = b['rgdp'].values[-1]
    r = {}
    tot = 0.0
    for k, s in enumerate(ALL6):
        r[(s, 'VA %')] = pct(m, 'va%d' % k)[-1]
        r[(s, 'L %')] = pct(m, 'l%d' % k)[-1]
        r[(s, 'K %')] = pct(m, 'k%d' % k)[-1]
        r[(s, 'VA/L %')] = pct(m, 'lp%d' % k)[-1]
        c = 100 * (x['va%d' % k].values[-1] - b['va%d' % k].values[-1]) / rg0
        r[(s, 'contrib pp')] = c
        tot += c
    r[('Total', 'real GDP %')] = pct(m, 'rgdp')[-1]
    r[('Total', 'residual pp')] = r[('Total', 'real GDP %')] - tot
    rows[name] = r
va = pd.DataFrame(rows).T
va.columns = pd.MultiIndex.from_tuples(va.columns)
# LaTeX: contributions and VA per worker only (the CSV has everything)
vt = pd.DataFrame({**{'%s' % s: va[(s, 'contrib pp')] for s in ALL6},
                   'Residual': va[('Total', 'residual pp')], 'Real GDP': va[('Total', 'real GDP %')]})
hdr = [' & \\multicolumn{6}{c}{Contribution of sectoral real VA to real GDP change (pp)} & & \\\\',
       '\\cmidrule(lr){2-7}',
       ' & ' + ' & '.join(ALL6) + ' & Residual & Real GDP (\\%) \\\\']
write(vt, 'tab_va_sectors_lr', '%.2f', hdr)
va.to_csv(out('tab_va_sectors_lr') + '.csv')   # overwrite write()'s CSV with the full version

# ---------------------------------------------------------------------------------------------
# Paths: labor decomposition and external accounts

lab, ext = [], []
for name, m in POLICIES.items():
    b, x = m
    T = len(x)
    dlg = 100 * np.log(x['lg'] / b['lg']).values
    dll = 100 * np.log(x['ll'] / b['ll']).values
    lab.append(pd.DataFrame({'policy': name, 'period': x['period'].values,
                             'dlog_goods': dlg, 'dlog_total': dll, 'dlog_goods_share': dlg - dll}))
    e = {'policy': name, 'period': x['period'].values}
    for lab_, mm in [('base', b), ('tariff', x)]:
        e['tb_y_' + lab_] = 100 * mm['nx_y'].values
        e['tb_goods_y_' + lab_] = 100 * mm['nxg_y'].values
        e['tb_svcs_y_' + lab_] = 100 * mm['nxsvc_y'].values
        e['tb_chn_y_' + lab_] = 100 * mm['nx1_y'].values
        e['tb_row_y_' + lab_] = 100 * mm['nx2_y'].values
        if 'b_y' in mm:
            e['nfa_y_' + lab_] = 100 * mm['b_y'].values
        if 'nc' in mm:
            e['c_y_' + lab_] = 100 * (mm['nc'] / mm['ngdp']).values
        e['i_y_' + lab_] = 100 * (mm['ni'] / mm['ngdp']).values
    e['real_exports_pct'] = pct(m, 'rex')
    e['real_imports_pct'] = pct(m, 'rim')
    e['real_goods_exports_pct'] = pct(m, 'rexg')
    e['real_goods_imports_pct'] = pct(m, 'rimg')
    e['real_svcs_exports_pct'] = pct(m, 'rexs4')
    e['real_svcs_imports_pct'] = pct(m, 'rims4')
    ext.append(pd.DataFrame(e))
    if 'nc' in x:
        # accounting check: TB = Y - C - I at current prices
        gap = (x['ngdp'] - x['nc'] - x['ni'] - x['nx']).abs().max() / x['ngdp'].mean()
        if gap > 1e-8:
            print('tables.py: WARNING %s: |Y - C - I - NX|/Y = %.2e' % (name, gap))
pd.concat(lab).to_csv(out('paths_labor_decomp') + '.csv', index=False)
pd.concat(ext).to_csv(out('paths_external') + '.csv', index=False)

print('tables.py: wrote tables to %s%s' % (args.out, ' (GHH suffix)' if args.ghh else ''))
