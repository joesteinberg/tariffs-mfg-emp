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
# load the results

def load_results(t,s,c,r,d,a):

    models_usa=[]
    models_usa.append(pd.read_csv('../c/output/vars0_usa_a%d.csv'%a))
    models_usa.append(pd.read_csv('../c/output/vars1_usa_t%d_s%d_c%d_r%d_d%d_a%d.csv'%(t,s,c,r,d,a)))
                        
    for m in models_usa:
        m['nx']=m['nx1']+m['nx2']
        m['lg']=m['l0']+m['l1']+m['l2']+m['l3']
        for s in [0,1,2,3,4]:
            m['nxs'+str(s)]=m['nxs'+str(s)+'-1']+m['nxs'+str(s)+'-2']
            m['exs'+str(s)]=m['exs'+str(s)+'-1']+m['exs'+str(s)+'-2']
            m['ims'+str(s)]=m['ims'+str(s)+'-1']+m['ims'+str(s)+'-2']
            
        m['nxsg'] = m['nxs0'] + m['nxs1'] + m['nxs2'] + m['nxs3']
        for col in m.columns:
            if('nx' in col or 'ex' in col or 'im' in col):
                m[col] = m[col]/m.ngdp                

    return models_usa


t=25
d=0
r=0
a=4

df_all = [load_results(t,s,2,r,d,a) for s in [0,1,3,4,6]]
df_chn = [load_results(t,s,0,r,d,a) for s in [0,1,3,4,6]]
df_row = [load_results(t,s,1,r,d,a) for s in [0,1,3,4,6]]
df_m = [load_results(t,s,2,r,d,a) for s in [10,11,13,14,16]]
df_f = [load_results(t,s,2,r,d,a) for s in [20,21,23,24,26]]

df_all_ret = load_results(t,6,2,1,0,4)
df_all_trans = load_results(t,6,2,0,1,4)
df_all_noadj = load_results(t,6,2,0,0,0)


slabels=['Oil','Steel','Toys','Cars']

def pct_chg(df,variable):
    return (100*(df[1][variable]/df[0][variable] - 1)).values

def diff_pct_lg(df,variable):
    return (100*(df[1][variable]-df[0][variable])/df[0]['lg']).values

def diff_pct_ll(df,variable):
    return (100*(df[1][variable]-df[0][variable])/df[0]['ll']).values

def diff_pct_ii(df,variable):
    return (100*(df[1][variable]-df[0][variable])/df[0]['ii']).values

def pp_chg(df,variable):
    return (100*(df[1][variable]-df[0][variable])).values

def create_fig(size=(3.25,3.4)):
    fig, ax = plt.subplots(figsize=size,tight_layout = {'pad': 0})
    ax.tick_params(axis='both', labelsize=10)
    ax.yaxis.label.set_size(10)
    sns.despine()
    return fig,ax


##########################################################################
# transitions: one set of plots per scenario

mperiods = df_all[0][0].period.values

def plot_dyn(data,fname,leg=1,ylim=None,lr_bars=False,ylabel=None,size=(3.25,3.4)):

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

    if(leg):
        ax.legend(frameon=True, loc='best',fontsize=8, framealpha=1, edgecolor='white',borderpad=0,borderaxespad=1)

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

    ax.set_xlabel('period')

    if ylabel !=None:
        ax.set_ylabel(ylabel)

    ax.axhline(0.0,color='black',lw=1)
            
    fig.tight_layout()
    plt.savefig('output/' + fname + '.pdf')
    plt.clf()

# sector-level employment in each scenario
for s in range(5):
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = diff_pct_lg(df_all[s],'l%d'%s2)
    data['Total'] = pct_chg(df_all[s],'lg')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_labor_goods_s%d'%s
    plot_dyn(data,fname,ylim=(-2,5.5),lr_bars=True,ylabel='change (percent total goods emp)',size=(3.25,3.4))        
    plt.close('all')

# go in each scenario
for s in range(5):
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pct_chg(df_all[s],'y%d'%s2)
    #data['Total'] = pct_chg(df_all[s],'lg')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_y_goods_s%d'%s
    plot_dyn(data,fname,lr_bars=True,ylabel='change in output (percent)',size=(3.25,3.4))
    plt.close('all')

# inv in each scenario
for s in range(5):
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pct_chg(df_all[s],'i%d'%s2)
    #data['Total'] = pct_chg(df_all[s],'ii')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_inv_goods_s%d'%s
    plot_dyn(data,fname,lr_bars=True,ylabel='change in investment (percent)',size=(3.25,3.4))
    plt.close('all')

# capital costs in each scenario
for s in range(5):
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pct_chg(df_all[s],'rks%d'%s2)
    #data['Total'] = pct_chg(df_all[s],'lg')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_rks_goods_s%d'%s
    plot_dyn(data,fname,lr_bars=True,ylabel='change in relative capital cost (percent)',size=(3.25,3.4),ylim=(-5,10))
    plt.close('all')

# materials costs in each scenario
for s in range(5):
    data = {}
    for s2 in range(4):
        data[slabels[s2]] = pct_chg(df_all[s],'pm2s%d'%s2)
    #data['Total'] = pct_chg(df_all[s],'lg')
    data = pd.DataFrame(data)
    fname = 'fig_dyn_pm2_goods_s%d'%s
    plot_dyn(data,fname,lr_bars=True,ylabel='change in relative material cost (percent)',size=(3.25,3.4),ylim=(-5,10))
    plt.close('all')
    
# GDP across scenarios
scenarios = ['Oil','Steel','Toys','Cars','All']
data={}
for s in [0,1,2,3,4]:
    data[scenarios[s]] = pct_chg(df_all[s],'rgdp')
fname = 'fig_dyn_macro'
data = pd.DataFrame(data)
plot_dyn(data,fname,ylabel = 'percent change',size=(3.25*0.5/0.4,3.4),lr_bars=True)

## resized main graph
data = {}
for s2 in range(4):
    data[slabels[s2]] = diff_pct_lg(df_all[4],'l%d'%s2)
data['Total'] = pct_chg(df_all[4],'lg')
data = pd.DataFrame(data)
fname = 'fig_dyn_labor_goods_atb_resize'
plot_dyn(data,fname,ylim=(-2,3),lr_bars=True,ylabel='change in sectoral emp (percent total goods emp)',size=(3.25,3.4))

# ATB tariffs w/ retaliation
data = {}
for s2 in range(4):
    data[slabels[s2]] = diff_pct_lg(df_all_ret,'l%d'%s2)
data['Total'] = pct_chg(df_all_ret,'lg')
data = pd.DataFrame(data)
fname = 'fig_dyn_labor_goods_atb_retaliation'
plot_dyn(data,fname,ylim=(-2,3),lr_bars=True,ylabel='change in sectoral emp (percent total goods emp)',size=(3.25,3.4))

# ATB tariffs w/ no adj costs
data = {}
for s2 in range(4):
    data[slabels[s2]] = diff_pct_lg(df_all_noadj,'l%d'%s2)
data['Total'] = pct_chg(df_all_noadj,'lg')
data = pd.DataFrame(data)
fname = 'fig_dyn_labor_goods_atb_noadj'
plot_dyn(data,fname,ylim=(-2,3),lr_bars=True,ylabel='change in sectoral emp (percent total goods emp)',size=(3.25,3.4))

# ATB tariffs w/ reversal
data = {}
for s2 in range(4):
    data[slabels[s2]] = diff_pct_lg(df_all_trans,'l%d'%s2)
data['Total'] = pct_chg(df_all_trans,'lg')
data = pd.DataFrame(data)
fname = 'fig_dyn_labor_goods_atb_temp'
plot_dyn(data,fname,ylim=(-2,3),lr_bars=True,ylabel='change in sectoral emp (percent total goods emp)',size=(3.25,3.4))


plt.close('all')
    


##########################################################################
# long-run aggregate goods-sector employment

width=0.75

data1 = np.array([pct_chg(df_all[s],'lg')[-1] for s in range(0,5)])
data_df1 = pd.DataFrame({'emp_chg1':data1,'scenario':scenarios}).sort_values(by='emp_chg1',ascending=False)

data2 = np.array([pct_chg(df_chn[s],'lg')[-1] for s in range(0,5)])
data_df2 = pd.DataFrame({'emp_chg2':data2,'scenario':scenarios})

data3 = np.array([pct_chg(df_row[s],'lg')[-1] for s in range(0,5)])
data_df3 = pd.DataFrame({'emp_chg3':data3,'scenario':scenarios})

data4 = np.array([pct_chg(df_all[s],'rgdp')[-1] for s in range(0,5)])
data_df4 = pd.DataFrame({'gdp_chg':data4,'scenario':scenarios})

data5 = np.array([pct_chg(df_m[s],'lg')[-1] for s in range(0,5)])
data_df5 = pd.DataFrame({'emp_chg_m':data5,'scenario':scenarios})

data6 = np.array([pct_chg(df_f[s],'lg')[-1] for s in range(0,5)])
data_df6 = pd.DataFrame({'emp_chg_f':data6,'scenario':scenarios})

data_df = pd.merge(left=data_df1,right=data_df2,how='left',on='scenario')
data_df = pd.merge(left=data_df,right=data_df3,how='left',on='scenario')
data_df = pd.merge(left=data_df,right=data_df4,how='left',on='scenario')
data_df = pd.merge(left=data_df,right=data_df5,how='left',on='scenario')
data_df = pd.merge(left=data_df,right=data_df6,how='left',on='scenario')

for s2 in range(0,4):
    data7 = np.array([diff_pct_lg(df_all[s],'l%d'%s2)[-1] for s in range(0,5)])
    data_df7 = pd.DataFrame({'emp_chg_s%d'%s2:data7,'scenario':scenarios})
    data_df = pd.merge(left=data_df,right=data_df7,how='left',on='scenario')

    
# just total employment
fig,ax=create_fig()
ax.bar(data_df.scenario, data_df.emp_chg1, width, color='black', alpha=0.99, edgecolor='black', lw=1)
ax.axhline(0.0,color='black',ls='-',lw=0.5)
ax.set_ylim(-2.5,5.5)
ax.set_xlabel('targeted sector')
ax.set_ylabel('percent chg in goods emp')
fig.tight_layout()    
plt.savefig('output/fig_lr_emp_across_scenarios_all.pdf')
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
ax.set_ylim(-2.5,5.5)
ax.set_xlabel('Targeted sector')
ax.set_ylabel('change (percent total goods emp)')
fig.tight_layout()    
plt.savefig('output/fig_lr_emp_across_scenarios_all_vs_sectors.pdf')
plt.clf()

# china only
fig,ax=create_fig()    
ax.bar(data_df.scenario, data_df.emp_chg2, width, color='black', alpha=0.99, edgecolor='black', lw=1)
ax.axhline(0.0,color='black',ls='-',lw=0.5)
ax.set_ylim(-2,3.4)
ax.set_xlabel('targeted sector')
ax.set_ylabel('percent chg in total goods emp')
fig.tight_layout()    
plt.savefig('output/fig_lr_emp_across_scenarios_chn.pdf')
plt.clf()

# row only
fig,ax=create_fig()    
ax.bar(data_df.scenario, data_df.emp_chg3, width, color='black', alpha=0.99, edgecolor='black', lw=1)
ax.axhline(0.0,color='black',ls='-',lw=0.5)
ax.set_ylim(-2,3.4)
ax.set_xlabel('targeted sector')
ax.set_ylabel('percent chg in total goods emp')
fig.tight_layout()    
plt.savefig('output/fig_lr_emp_across_scenarios_row.pdf')
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
ax.set_ylim(-2,3.4)
ax.set_xlabel('targeted sector')
ax.set_ylabel('percent change in total goods emp')
fig.tight_layout()    
plt.savefig('output/fig_lr_emp_across_scenarios_compare.pdf')
plt.clf()

# compare emp vs gdp
fig,ax=create_fig()    
ax.bar(x-width/2, data_df.emp_chg1, width, color='black', alpha=0.99, edgecolor='black', lw=1, label='Goods employment')
ax.bar(x+width/2, data_df.gdp_chg, width, color=colors[0], alpha=0.99, edgecolor='black', lw=1, label='Real GDP')
ax.set_xticks(x, data_df.scenario)
ax.axhline(0.0,color='black',ls='-',lw=1)
ax.legend(frameon=True, loc='upper right',fontsize=8, framealpha=1, edgecolor='white',borderpad=0,borderaxespad=1)
ax.set_ylim(-2,3.4)
ax.set_xlabel('targeted sector')
ax.set_ylabel('percent change')
fig.tight_layout()    
plt.savefig('output/fig_lr_emp_vs_gdp.pdf')
plt.clf()


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
ax.set_ylim(-2,3.4)
ax.set_xlabel('targeted sector')
ax.set_ylabel('pct change in total goods emp')
fig.tight_layout()    
plt.savefig('output/fig_lr_emp_all_vs_m_vs_f.pdf')
plt.clf()



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
    ax.set_ylim(-1.5,5.5)
    fig.tight_layout()    
    plt.savefig('output/'+fname+'.pdf')
    plt.clf()
    

for s in range(5):
    fname='fig_lr_goods_labor_s%d'%s
    plot_lr_emp(df_all[s],fname)
    
plt.close('all')
