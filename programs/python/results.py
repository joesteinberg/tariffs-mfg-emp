import numpy as np
import pandas as pd
import matplotlib as mpl
import pandas as pd

import locale
locale.setlocale(locale.LC_ALL,'en_US.utf8')

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Roboto'] + plt.rcParams['font.sans-serif']
plt.rcParams['font.weight'] = 365

import seaborn as sns

colors=['#e41a1c','#377eb8','#984ea3','#4daf4a','black']
hatches=[None,'x','/','.','+']
lses = [(5,(10,3)),'dashdot',(0,(5,1)),(0,(3,1,1,1,1,1)),'solid']

def autolabel(ax,rects,sign):
    for rect in rects:
        height=sign*rect.get_height()
        ax.annotate('%0.2f' % float(height),
                    (rect.get_x()+rect.get_width()/2.0,height),
                    ha='center',
                    va='center',
                    xytext=(0,sign*20),
                    textcoords='offset points')

vmax = lambda v: [max(x,0.0) for x in v]
vmin = lambda v: [min(x,0.0) for x in v]

def which_bottom(bottom_p, bottom_n, data):
    bottom=np.zeros(len(data))
    for i in range(len(data)):
        bottom[i]=bottom_p[i] if data[i]>0.0 else bottom_n[i]
    return bottom


##########################################################################
# load the results (see model_io.py for options: --cdir, --out, --ghh, --pz, --partial)

import os
from model_io import (parse_args, load_results, load_results_trump2025, have, check_inputs,
                      note_skipped, pct_chg, diff_pct_lg, pp_chg, SLABELS, SCENARIOS, SCEN_CODES)

args = parse_args()

def savefig(fname):
    plt.savefig(os.path.join(args.out, fname + args.suffix + '.pdf'))

def skipped(what):
    # only reached with --partial or for figures outside the GHH set (check_inputs below)
    why = 'not part of the sensitivity run set' if args.sens else 'model output not in %s' % args.cdir
    note_skipped(args, 'results.py', '%s (%s)' % (what, why))

t=25
d=0
r=0
a=4

df_all = [load_results(t,s,2,r,d,a,args) for s in SCEN_CODES]
df_chn = [load_results(t,s,0,r,d,a,args) for s in SCEN_CODES]
df_row = [load_results(t,s,1,r,d,a,args) for s in SCEN_CODES]
df_m = [load_results(t,s,2,r,d,a,args) for s in [10,11,13,14,16]]
df_f = [load_results(t,s,2,r,d,a,args) for s in [20,21,23,24,26]]

df_all_ret = load_results(t,6,2,1,0,4,args)
df_all_trans = load_results(t,6,2,0,1,4,args)
df_all_noadj = load_results(t,6,2,0,0,0,args)
df_trump2025 = load_results_trump2025(0,0,4,args)
check_inputs(args, 'results.py')

slabels = SLABELS

def create_fig(size=(3.25,3.25)):
    fig, ax = plt.subplots(figsize=size,tight_layout = {'pad': 0})
    ax.tick_params(axis='both', labelsize=10)
    ax.yaxis.label.set_size(10)
    sns.despine()
    return fig,ax


##########################################################################
# transitions: one set of plots per scenario

_avail = [x for x in df_all + [df_trump2025] if x is not None]
if not _avail:
    raise SystemExit('results.py: no model output found in %s' % args.cdir)
mperiods = _avail[0][0].period.values

def _legend_clear(fig, ax, leg):
    """True if the legend box overlaps no labeled line and no bar."""
    from matplotlib.patches import Rectangle
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = leg.get_window_extent(r)
    for l in ax.get_lines():
        if l.get_label().startswith('_'):
            continue
        xy = ax.transData.transform(np.column_stack([l.get_xdata(), l.get_ydata()]))
        if len(xy) > 1:   # densify so segments crossing the box are caught
            t = np.linspace(0, 1, 10)[:, None]
            xy = np.vstack([xy[i] + (xy[i+1]-xy[i])*t for i in range(len(xy)-1)])
        if bb.count_contains(xy) > 0:
            return False
    for p in ax.patches:
        if isinstance(p, Rectangle) and p.get_height() != 0 and p.get_window_extent(r).overlaps(bb):
            return False
    return True

def place_legend(fig, ax, ncl):
    """Put the legend where it covers no data: try standard locations (with the requested and
    a wider column count), and if none is clear, add headroom above the data and put it there."""
    kw = dict(frameon=True, fontsize=8, framealpha=1, edgecolor='white', borderpad=0, borderaxespad=1)
    n = len([l for l in ax.get_lines() if not l.get_label().startswith('_')])
    for nc in dict.fromkeys([ncl, max(ncl, 2), max(ncl, 3)]):
        if nc > n:
            continue
        for loc in ['upper left', 'upper center', 'upper right', 'lower left', 'lower center',
                    'lower right', 'center left', 'center right', 'center']:
            leg = ax.legend(loc=loc, ncol=nc, **kw)
            if _legend_clear(fig, ax, leg):
                return
    leg = ax.legend(loc='upper left', ncol=max(ncl, 2), **kw)
    lo, hi = ax.get_ylim()
    for _ in range(20):
        if _legend_clear(fig, ax, leg):
            return
        hi = hi + 0.1*(hi - lo)
        ax.set_ylim(lo, hi)

def plot_dyn(data,fname,leg=1,ylim=None,lr_bars=False,ylabel=None,size=(3.25,3.25),noxlab=False,ncl=1):

    maxT=51
    W=2
    
    fig,ax=create_fig(size)

    n = len(data.columns)

    lns=[]
    cnt=0
    for col in data.columns:
        l = ax.plot(mperiods[0:maxT], data[col][0:maxT].values,
                    linestyle=lses[cnt], color=colors[cnt],
                    marker=None,alpha=0.99,zorder=0,linewidth=1,label=col)
        lns.append(l)
            
        cnt=cnt+1

    if(lr_bars):
        ax.bar(maxT+W*n/2 + 2*np.arange(n), [data[col].values[-1] for col in data.columns],
               width=W, alpha=0.99, edgecolor='black', lw=1, color=colors)
        locs, labels = plt.xticks()
        ax.set_xticks([0,10,20,30,40,50,60],labels=['0','10','20','30','40','50','LR'])
        #xticks = ax.get_xticks()
        #xlabels = [item.get_text() for item in ax.get_xticklabels()]
        #xlabels[-1] = 'LR'
        #ax.set_xticks(xticks, labels=xlabels)
        
    if ylim !=None:
        ax.set_ylim(ylim)

    if(noxlab==False):
        ax.set_xlabel('period')

    if ylabel !=None:
        ax.set_ylabel(ylabel)

    ax.axhline(0.0,color='black',lw=1)
    fig.tight_layout()

    if(leg):
        place_legend(fig, ax, ncl)

    savefig(fname)
    plt.clf()

# industry-level goods employment in each scenario
for s in [k for k in range(5) if df_all[k] is not None]:
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = diff_pct_lg(df_all[s],'l%d'%s2)
    data['Total'] = pct_chg(df_all[s],'lg')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_labor_goods_s%d'%s
    plot_dyn(data,fname,ylim=(-2.5,5.5),lr_bars=True,ylabel='change (percent total goods emp)',size=(3.25,3.25))        
    plt.close('all')

# sector-level employment in each scenario
for s in [k for k in range(5) if df_all[k] is not None]:
    data = {}
    data['Goods'] = pct_chg(df_all[s],'lg')
    data['Svcs'] = pct_chg(df_all[s],'l4')
    data['Const'] = pct_chg(df_all[s],'l5')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_labor_sectors_s%d'%s
    plot_dyn(data,fname,leg=(s==2),ylim=(-3,3.5),lr_bars=True,size=(2.5,1.8),noxlab=True)
    plt.close('all')
    
# go in each scenario
for s in [k for k in range(5) if df_all[k] is not None]:
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pct_chg(df_all[s],'y%d'%s2)
    #data['Total'] = pct_chg(df_all[s],'lg')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_y_goods_s%d'%s
    plot_dyn(data,fname,leg=(s==2),lr_bars=True,ylabel=None,size=(2.5,1.8))
    plt.close('all')

# inv in each scenario
for s in [k for k in range(5) if df_all[k] is not None]:
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pct_chg(df_all[s],'i%d'%s2)
    #data['Total'] = pct_chg(df_all[s],'ii')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_inv_goods_s%d'%s
    plot_dyn(data,fname,lr_bars=True,ylabel='change in investment (percent)',size=(3.25,3.25))
    plt.close('all')

# capital costs in each scenario
for s in [k for k in range(5) if df_all[k] is not None]:
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pct_chg(df_all[s],'rks%d'%s2)
    #data['Total'] = pct_chg(df_all[s],'lg')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_rks_goods_s%d'%s
    plot_dyn(data,fname,lr_bars=True,ylabel='change in relative capital cost (percent)',size=(3.25,3.25),ylim=(-10,6))
    plt.close('all')

# materials costs in each scenario
for s in [k for k in range(5) if df_all[k] is not None]:
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pct_chg(df_all[s],'pm2s%d'%s2)
    #data['Total'] = pct_chg(df_all[s],'lg')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_pm2_goods_s%d'%s
    plot_dyn(data,fname,lr_bars=True,ylabel='change in relative material cost (percent)',size=(3.25,3.25),ylim=(-1,5))
    plt.close('all')

# intermediate prices in each scenario
for s in [k for k in range(5) if df_all[k] is not None]:
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pct_chg(df_all[s],'pm%d'%s2)
    #data['Total'] = pct_chg(df_all[s],'lg')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_pm_goods_s%d'%s
    #plot_dyn(data,fname,ylim=(-2,4),lr_bars=True,size=(2.5,1.8),noxlab=True)
    plot_dyn(data,fname,leg=(s==2),lr_bars=True,ylabel=None,size=(2.5,1.8),ylim=(-5,20),ncl=2)
    plt.close('all')

# final prices in each scenario
for s in [k for k in range(5) if df_all[k] is not None]:
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pct_chg(df_all[s],'pf%d'%s2)
    #data['Total'] = pct_chg(df_all[s],'lg')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_pf_goods_s%d'%s
    plot_dyn(data,fname,leg=(s==2),lr_bars=True,ylabel=None,size=(2.5,1.8),ylim=(-3,10),ncl=2)
    plt.close('all')

# imports in each scenario
for s in [k for k in range(5) if df_all[k] is not None]:
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pct_chg(df_all[s],'ims%d'%s2)
    data['Total'] = pct_chg(df_all[s],'img')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_imports_goods_s%d'%s
    plot_dyn(data,fname,leg=(s==2),lr_bars=True,ylabel=None,size=(2.5,1.8),ylim=(-80,25),ncl=2)
    plt.close('all')

# exports in each scenario
for s in [k for k in range(5) if df_all[k] is not None]:
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pct_chg(df_all[s],'exs%d'%s2)
    data['Total'] = pct_chg(df_all[s],'exg')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_exports_goods_s%d'%s
    plot_dyn(data,fname,leg=(s==2),lr_bars=True,ylabel=None,size=(2.5,1.8),ylim=(-80,5),ncl=2)
    plt.close('all')

# net exports in each scenario
for s in [k for k in range(5) if df_all[k] is not None]:
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pp_chg(df_all[s],'nxs%d_y'%s2)
    data['Total'] = pp_chg(df_all[s],'nxg_y')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_nx_goods_s%d'%s
    plot_dyn(data,fname,leg=(s==2),lr_bars=True,ylabel=None,size=(2.5,1.8),ylim=(-0.75,2),ncl=2)
    plt.close('all')

# real exchange rate
for s in [k for k in range(5) if df_all[k] is not None]:
    data = {}
    data['REER'] = pct_chg(df_all[s],'rer')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_rer_s%d'%s
    plot_dyn(data,fname,leg=False,lr_bars=True,ylabel=None,size=(2.5,1.8),ylim=(-20,2.5))
    plt.close('all')


    
    
# GDP and consumption across scenarios
scenarios = SCENARIOS
if have(df_all):
    for var,fname in [('rgdp','fig_dyn_macro'),('cagg','fig_dyn_consumption')]:
        data = pd.DataFrame({scenarios[s]: pct_chg(df_all[s],var) for s in range(5)})
        plot_dyn(data,fname,ylabel = 'percent change',size=(3.25*0.5/0.4,3.4),lr_bars=True)
else:
    skipped('fig_dyn_macro, fig_dyn_consumption')

def plot_goods_labor(df, fname, ylabel='change (percent total goods emp)'):
    if df is None:
        skipped(fname)
        return
    data = {slabels[s2]: diff_pct_lg(df,'l%d'%s2) for s2 in range(4)}
    data['Total'] = pct_chg(df,'lg')
    plot_dyn(pd.DataFrame(data),fname,ylim=(-2.5,2),lr_bars=True,ylabel=ylabel,size=(3.25,3.25))

plot_goods_labor(df_all[4], 'fig_dyn_labor_goods_atb_resize', 'change in sectoral emp (percent total goods emp)')
plot_goods_labor(df_all_ret, 'fig_dyn_labor_goods_atb_retaliation')   # ATB tariffs w/ retaliation
plot_goods_labor(df_all_noadj, 'fig_dyn_labor_goods_atb_noadj')       # ATB tariffs w/ no adj costs
plot_goods_labor(df_all_trans, 'fig_dyn_labor_goods_atb_temp')        # ATB tariffs w/ reversal
plot_goods_labor(df_trump2025, 'fig_dyn_labor_goods_trump2025')       # Trump 2025 tariff matrix

plt.close('all')
    


##########################################################################
# external accounts and labor decomposition

def plot_tb(df, fname):
    """Change in the trade balance / GDP (pp): by type, bilateral, and total."""
    if df is None:
        skipped(fname)
        return
    data = pd.DataFrame({'Goods': pp_chg(df,'nxg_y'), 'Services': pp_chg(df,'nxsvc_y'),
                         'China': pp_chg(df,'nx1_y'), 'RoW': pp_chg(df,'nx2_y'), 'Total': pp_chg(df,'nx_y')})
    plot_dyn(data,fname,lr_bars=True,ylabel='change in trade balance (pp of GDP)',size=(3.25,3.25),ncl=2)
    plt.close('all')

for s in range(5):
    plot_tb(df_all[s], 'fig_dyn_tb_s%d'%s)
plot_tb(df_trump2025, 'fig_dyn_tb_trump2025')

if have(df_all):
    for var,fname,ylabel in [('nx_y','fig_dyn_tb_across','change in trade balance (pp of GDP)'),
                             ('b_y','fig_dyn_nfa_across','change in NFA (pp of GDP)')]:
        data = pd.DataFrame({scenarios[s]: pp_chg(df_all[s],var) for s in range(5)})
        plot_dyn(data,fname,ylabel=ylabel,size=(3.25*0.5/0.4,3.4),lr_bars=True)
        plt.close('all')
else:
    skipped('fig_dyn_tb_across, fig_dyn_nfa_across')

def plot_labor_decomp(df, fname):
    """dlog goods emp = dlog total labor + dlog goods share (log points x 100)."""
    if df is None:
        skipped(fname)
        return
    g = 100*np.log(df[1]['lg']/df[0]['lg']).values
    l = 100*np.log(df[1]['ll']/df[0]['ll']).values
    data = pd.DataFrame({'Total labor': l, 'Goods share': g-l, 'Goods emp': g})
    plot_dyn(data,fname,lr_bars=True,ylabel='log change (x100)',size=(3.25,3.25))
    plt.close('all')

plot_labor_decomp(df_all[4], 'fig_dyn_labor_decomp_s4')
plot_labor_decomp(df_trump2025, 'fig_dyn_labor_decomp_trump2025')


##########################################################################
# long-run aggregate goods-sector employment

width=0.75

lr = lambda dfs, var, f=pct_chg: np.array([f(x,var)[-1] if x is not None else np.nan for x in dfs])

data_df = pd.DataFrame({'scenario':scenarios,
                        'emp_chg1':lr(df_all,'lg'), 'emp_chg2':lr(df_chn,'lg'), 'emp_chg3':lr(df_row,'lg'),
                        'gdp_chg':lr(df_all,'rgdp'), 'emp_chg_m':lr(df_m,'lg'), 'emp_chg_f':lr(df_f,'lg')})
for s2 in range(0,4):
    data_df['emp_chg_s%d'%s2] = lr(df_all,'l%d'%s2,diff_pct_lg)
data_df = data_df.sort_values(by='emp_chg1',ascending=False).reset_index(drop=True)

    
if have(df_all):
    # just total employment
    fig,ax=create_fig()
    ax.bar(data_df.scenario, data_df.emp_chg1, width, color='black', alpha=0.99, edgecolor='black', lw=1)
    ax.axhline(0.0,color='black',ls='-',lw=0.5)
    ax.set_ylim(-3,6)
    ax.set_xlabel('targeted sector')
    ax.set_ylabel('percent chg in goods emp')
    fig.tight_layout()    
    savefig('fig_lr_emp_across_scenarios_all')
    plt.clf()

    # add stacked bar showing sectoral impacts
    x = np.arange(len(scenarios))
    width=0.75/2

    fig,ax=create_fig()    
    ax.bar(x-width/2, data_df.emp_chg1, width, color='black', alpha=0.99, edgecolor='black', lw=1, label='Total')

    bottom_pos=np.zeros(len(scenarios))
    bottom_neg=np.zeros(len(scenarios))

    bcnt=0
    for s2 in range(0,4):
        values = data_df['emp_chg_s%d'%s2]
        bottom=which_bottom(bottom_pos,bottom_neg,values)
        ax.bar(x+width/2, values, width, label=slabels[s2], bottom=bottom, color=colors[bcnt], alpha=0.99, edgecolor='black', lw=0.5)
        bottom_pos=bottom_pos+vmax(values)
        bottom_neg=bottom_neg+vmin(values)
        bcnt += 1

    ax.set_xticks(x, data_df.scenario)
    ax.axhline(0.0,color='black',ls='-',lw=1)
    ax.legend(frameon=True, loc='upper right',fontsize=8, framealpha=1, edgecolor='white',borderpad=0,borderaxespad=1)
    ax.set_ylim(-3,6)
    ax.set_xlabel('targeted sector')
    ax.set_ylabel('change (percent total goods emp)')
    fig.tight_layout()    
    savefig('fig_lr_emp_across_scenarios_all_vs_sectors')
    plt.clf()
else:
    skipped('LR employment bars across scenarios')

if have(df_all, df_chn, df_row):
    # china only
    fig,ax=create_fig()    
    ax.bar(data_df.scenario, data_df.emp_chg2, width, color='black', alpha=0.99, edgecolor='black', lw=1)
    ax.axhline(0.0,color='black',ls='-',lw=0.5)
    ax.set_ylim(-2.5,3.5)
    ax.set_xlabel('targeted sector')
    ax.set_ylabel('percent chg in total goods emp')
    fig.tight_layout()    
    savefig('fig_lr_emp_across_scenarios_chn')
    plt.clf()

    # row only
    fig,ax=create_fig()    
    ax.bar(data_df.scenario, data_df.emp_chg3, width, color='black', alpha=0.99, edgecolor='black', lw=1)
    ax.axhline(0.0,color='black',ls='-',lw=0.5)
    ax.set_ylim(-2.5,3.5)
    ax.set_xlabel('targeted sector')
    ax.set_ylabel('percent chg in total goods emp')
    fig.tight_layout()    
    savefig('fig_lr_emp_across_scenarios_row')
    plt.clf()

    # compare all vs china/row only
    x = np.arange(len(scenarios))
    width=0.75/2

    fig,ax=create_fig()    
    ax.bar(x-width/2, data_df.emp_chg1, width, color='black', alpha=0.99, edgecolor='black', lw=1, label='Both')

    bottom_pos=np.zeros(len(scenarios))
    bottom_neg=np.zeros(len(scenarios))

    bcnt=0
    values = data_df.emp_chg2
    bottom=which_bottom(bottom_pos,bottom_neg,values)
    ax.bar(x+width/2, values, width, label='China only', bottom=bottom, color=colors[bcnt], alpha=0.99, edgecolor='black', lw=0.5)
    bottom_pos=bottom_pos+vmax(values)
    bottom_neg=bottom_neg+vmin(values)
    bcnt += 1

    values = data_df.emp_chg3
    bottom=which_bottom(bottom_pos,bottom_neg,values)
    ax.bar(x+width/2, values, width, label='RoW only', bottom=bottom, color=colors[bcnt], alpha=0.99, edgecolor='black', lw=0.5)
    bottom_pos=bottom_pos+vmax(values)
    bottom_neg=bottom_neg+vmin(values)
    bcnt += 1

    ax.set_xticks(x, data_df.scenario)
    ax.axhline(0.0,color='black',ls='-',lw=1)
    ax.legend(frameon=True, loc='upper right',fontsize=8, framealpha=1, edgecolor='white',borderpad=0,borderaxespad=1)
    ax.set_ylim(-2.5,3.5)
    ax.set_xlabel('targeted sector')
    ax.set_ylabel('percent change in total goods emp')
    fig.tight_layout()    
    savefig('fig_lr_emp_across_scenarios_compare')
    plt.clf()
else:
    skipped('China-only / RoW-only comparison')

if have(df_all):
    # compare emp vs gdp
    fig,ax=create_fig()    
    ax.bar(x-width/2, data_df.emp_chg1, width, color='black', alpha=0.99, edgecolor='black', lw=1, label='Goods employment')
    ax.bar(x+width/2, data_df.gdp_chg, width, color=colors[0], alpha=0.99, edgecolor='black', lw=1, label='Real GDP')
    ax.set_xticks(x, data_df.scenario)
    ax.axhline(0.0,color='black',ls='-',lw=1)
    ax.legend(frameon=True, loc='upper right',fontsize=8, framealpha=1, edgecolor='white',borderpad=0,borderaxespad=1)
    ax.set_ylim(-2.5,3.5)
    ax.set_xlabel('targeted sector')
    ax.set_ylabel('percent change')
    fig.tight_layout()    
    savefig('fig_lr_emp_vs_gdp')
    plt.clf()
else:
    skipped('LR employment vs GDP')

if have(df_all, df_m, df_f):
    # add stacked m vs f
    x = np.arange(len(scenarios))
    width=0.75/2

    fig,ax=create_fig()    
    ax.bar(x-width/2, data_df.emp_chg1, width, color='black', alpha=0.99, edgecolor='black', lw=1, label='Both')

    bottom_pos=np.zeros(len(scenarios))
    bottom_neg=np.zeros(len(scenarios))

    bcnt=0
    values = data_df.emp_chg_m
    bottom=which_bottom(bottom_pos,bottom_neg,values)
    ax.bar(x+width/2, values, width, label='Intermediates only', bottom=bottom, color=colors[bcnt], alpha=0.99, edgecolor='black', lw=0.5)
    bottom_pos=bottom_pos+vmax(values)
    bottom_neg=bottom_neg+vmin(values)
    bcnt += 1

    values = data_df.emp_chg_f
    bottom=which_bottom(bottom_pos,bottom_neg,values)
    ax.bar(x+width/2, values, width, label='Final goods only', bottom=bottom, color=colors[bcnt], alpha=0.99, edgecolor='black', lw=0.5)
    bottom_pos=bottom_pos+vmax(values)
    bottom_neg=bottom_neg+vmin(values)
    bcnt += 1

    ax.set_xticks(x, data_df.scenario)
    ax.axhline(0.0,color='black',ls='-',lw=1)
    ax.legend(frameon=True, loc='upper right',fontsize=8, framealpha=1, edgecolor='white',borderpad=0,borderaxespad=1)
    ax.set_ylim(-2.5,3.5)
    ax.set_xlabel('targeted sector')
    ax.set_ylabel('pct change in total goods emp')
    fig.tight_layout()    
    savefig('fig_lr_emp_all_vs_m_vs_f')
    plt.clf()
else:
    skipped('intermediates-only / final-only comparison')



##########################################################################
# long-run goods-sector employment breakdown

width=0.75

def plot_lr_emp(df,fname):

    fig,ax=create_fig()
    ax.bar(slabels, [diff_pct_lg(df,'l%d'%s2)[-1] for s2 in range(4)],
                     width, alpha=0.99, edgecolor='black', lw=1, color=colors)

    ax.axhline(0.0,color='black',ls='-',lw=0.5)
    ax.set_xlabel('Sector')
    ax.set_ylabel('change in sectoral emp (percent total goods emp)')
    ax.set_ylim(-2.5,5.5)
    fig.tight_layout()    
    savefig(fname)
    plt.clf()
    

for s in [k for k in range(5) if df_all[k] is not None]:
    fname='fig_lr_goods_labor_s%d'%s
    plot_lr_emp(df_all[s],fname)
    
plt.close('all')
