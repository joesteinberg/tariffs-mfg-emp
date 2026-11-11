"""Figure 3(c): long-run goods employment under the baseline, Caliendo--Parro, and KRS Armington
elasticities, for each single-sector tariff and tariffs on all goods.

usage: python3 fig_sig.py [--out output]
Reads output/tab_sig.csv (written by sig_table.py); writes fig_lr_emp_sig.pdf. Styling follows
the long-run bar charts in results.py; bars are ordered by the baseline effect, as there.
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Roboto'] + plt.rcParams['font.sans-serif']
plt.rcParams['font.weight'] = 365

ap = argparse.ArgumentParser()
ap.add_argument('--out', default='output')
args = ap.parse_args()

tab = pd.read_csv(os.path.join(args.out, 'tab_sig.csv')).set_index('vector')
CASES = [('base', 'Baseline', 'black'), ('cp', 'Caliendo-Parro', '#377eb8'), ('krs', 'KRS', '#e41a1c')]
POL = sorted(['Oil', 'Steel', 'Toys', 'Cars', 'All'], key=lambda p: -tab.loc['base', 'emp_%s' % p])

x = np.arange(len(POL))
w = 0.75 / len(CASES)
fig, ax = plt.subplots(figsize=(3.25, 3.25), tight_layout={'pad': 0})
ax.tick_params(axis='both', labelsize=10)
ax.yaxis.label.set_size(10)
sns.despine()
for k, (v, lab, col) in enumerate(CASES):
    ax.bar(x + (k - (len(CASES) - 1) / 2) * w, [tab.loc[v, 'emp_%s' % p] for p in POL], w,
           color=col, alpha=0.99, edgecolor='black', lw=0.5, label=lab)
ax.set_xticks(x, POL)
ax.axhline(0.0, color='black', ls='-', lw=1)
vals = np.array([[tab.loc[v, 'emp_%s' % p] for p in POL] for v, _, _ in CASES])
ax.set_ylim(min(-3.5, vals.min() - 0.5), max(8, vals.max() + 0.6))
ax.legend(frameon=True, loc='upper right', fontsize=8, framealpha=1, edgecolor='white', borderpad=0, borderaxespad=1)
ax.set_xlabel('targeted sector')
ax.set_ylabel('pct change in total goods emp')
fig.tight_layout()
plt.savefig(os.path.join(args.out, 'fig_lr_emp_sig.pdf'))
