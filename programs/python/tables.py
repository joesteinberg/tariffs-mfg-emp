"""Tables and decompositions from the model output (CSV + LaTeX tabular snippets).

usage: python3 tables.py [--cdir ../c/output] [--out output] [--ghh] [--partial]

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
from model_io import (parse_args, load_results, load_results_trump2025, check_inputs, note_skipped,
                      SLABELS, SCEN_CODES)

args = parse_args()

POLICIES = {'Oil': None, 'Steel': None, 'Toys': None, 'Cars': None, 'All': None, '2025': None}
for name, s in zip(['Oil', 'Steel', 'Toys', 'Cars', 'All'], SCEN_CODES):
    POLICIES[name] = load_results(25, s, 2, 0, 0, 4, args)
POLICIES['2025'] = load_results_trump2025(0, 0, 4, args)
check_inputs(args, 'tables.py')
missing = [k for k, v in POLICIES.items() if v is None]
if missing:
    note_skipped(args, 'tables.py', 'rows for %s (model output not in %s)' % (', '.join(missing), args.cdir))
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
        if args.nc == 3:
            e['tb_row_y_' + lab_] = 100 * mm['nx2_y'].values
        else:   # rest of world = low- + high-tariff partners
            e['tb_row_y_' + lab_] = 100 * (mm['nx2_y'] + mm['nx3_y']).values
            e['tb_ltp_y_' + lab_] = 100 * mm['nx2_y'].values
            e['tb_htp_y_' + lab_] = 100 * mm['nx3_y'].values
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

# ---------------------------------------------------------------------------------------------
# 4 regions: trade diversion. For each targeted sector (and all goods) and each set of targeted
# partners (China, low-tariff, high-tariff partners only; all partners): LR goods employment (%)
# and LR change in US real imports of the targeted goods from each region (%).

if args.nc == 4:
    TGT = {'Oil': [0], 'Steel': [1], 'Toys': [2], 'Cars': [3], 'All': [0, 1, 2, 3]}
    CASES = [('China', 0), ('Low-tariff', 3), ('High-tariff', 4), ('All', 2)]
    rows = []
    for name, s in zip(['Oil', 'Steel', 'Toys', 'Cars', 'All'], SCEN_CODES):
        for lab, c in CASES:
            m = load_results(25, s, c, 0, 0, 4, args)
            if m is None:
                note_skipped(args, 'tables.py', 'tab_diversion row %s / %s' % (name, lab))
                continue
            b, x = m
            imp = lambda d, j: sum(d['rims%d-%d' % (k, j)].values[-1] for k in TGT[name])
            r = {'targeted sector': name, 'tariffs on': lab, 'goods emp': pct(m, 'lg')[-1]}
            for j, reg in [(1, 'China'), (2, 'Low-tariff'), (3, 'High-tariff')]:
                r['imports from ' + reg] = 100 * (imp(x, j) / imp(b, j) - 1)
            rows.append(r)
    # the regions: share of US imports of each goods sector (no-tariff equilibrium = 2020 ICIO)
    # and 2025 tariff increase by sector (the 2025 run's tariff schedule)
    if '2025' in POLICIES:
        b0, x25 = POLICIES['2025']
        bl = []
        for j, reg in [(1, 'China'), (2, 'Low-tariff partners'), (3, 'High-tariff partners')]:
            r = {'region': reg}
            for k, sec in enumerate(SLABELS):
                tot = sum(b0['ims%d-%d' % (k, jj)].values[0] for jj in (1, 2, 3))
                r['import share ' + sec] = 100 * b0['ims%d-%d' % (k, j)].values[0] / tot
            totg = sum(b0['ims%d-%d' % (k, jj)].values[0] for k in range(4) for jj in (1, 2, 3))
            r['import share goods'] = 100 * sum(b0['ims%d-%d' % (k, j)].values[0] for k in range(4)) / totg
            for k, sec in enumerate(SLABELS):
                r['tariff ' + sec] = 100 * x25['taus%d-%d' % (k, j)].values[-1]
            # import-weighted (no-tariff imports) average over the four goods sectors
            w = np.array([b0['ims%d-%d' % (k, j)].values[0] for k in range(4)])
            r['tariff goods'] = 100 * np.dot(w, [x25['taus%d-%d' % (k, j)].values[-1] for k in range(4)]) / w.sum()
            bl.append(r)
        bl = pd.DataFrame(bl)
        bl.to_csv(out('tab_blocs') + '.csv', index=False)
        secs = SLABELS + ['goods']
        lines = ['\\begin{tabular}{l' + 'r' * 10 + '}', '\\toprule',
                 ' & \\multicolumn{5}{c}{Share of U.S. imports (\\%)} & \\multicolumn{5}{c}{2025 tariff increase (\\%)} \\\\',
                 '\\cmidrule(lr){2-6}\\cmidrule(lr){7-11}',
                 ' & ' + ' & '.join(['``%s\'\'' % s if s != 'goods' else 'All goods' for s in secs] * 2) + ' \\\\',
                 '\\midrule']
        for _, r in bl.iterrows():
            lines.append(r['region'] + ' & ' + ' & '.join('%.1f' % r['import share ' + s] for s in secs)
                         + ' & ' + ' & '.join('%.1f' % r['tariff ' + s] for s in secs) + ' \\\\')
        lines += ['\\bottomrule', '\\end{tabular}']
        with open(out('tab_blocs') + '.tex', 'w') as fh:
            fh.write('\n'.join(lines) + '\n')

    if rows:
        dv = pd.DataFrame(rows)
        dv.to_csv(out('tab_diversion') + '.csv', index=False)
        lines = ['\\begin{tabular}{llrrrr}', '\\toprule',
                 ' & & Goods & \\multicolumn{3}{c}{Real imports of targeted goods from} \\\\',
                 '\\cmidrule(lr){4-6}',
                 'Targeted sector & Tariffs on & employment & China & Low-tariff & High-tariff \\\\']
        for i, r in dv.iterrows():
            if i == 0 or r['targeted sector'] != dv['targeted sector'][i - 1]:
                lines.append('\\midrule')
                first = r['targeted sector']
            else:
                first = ''
            lines.append('%s & %s & %.2f & %.1f & %.1f & %.1f \\\\' % (
                first, r['tariffs on'], r['goods emp'], r['imports from China'],
                r['imports from Low-tariff'], r['imports from High-tariff']))
        lines += ['\\bottomrule', '\\end{tabular}']
        with open(out('tab_diversion') + '.tex', 'w') as fh:
            fh.write('\n'.join(lines) + '\n')

print('tables.py: wrote tables to %s%s' % (args.out, (' (suffix %s)' % args.suffix) if args.suffix else ''))
