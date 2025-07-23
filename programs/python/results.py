import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':16})
mpl.rc('font',size=16)
mpl.rc('text', usetex=True)
mpl.rc('lines',linewidth=1.5)
mpl.rc('savefig',bbox='tight')
mpl.rc('savefig',format='pdf')
mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=10)

hatches=[None,'x','/','.','+']
colors=['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
dashes = [(None,None),(10,2),(3,1)]
#colors=['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854']
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

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

tau=25
c=0
d=0
r=0

#sprintf(fname2,"output/%s_t%d_s%d_c%d_r%d_d%d_a%d.csv",
#	fname,(int)(tariff*100),target_sector_flag,target_country_flag,retaliation_flag,duration_flag,adjustment_flag);


def load_results(a):
    
    models_usa = []
    models_usa.append(pd.read_csv('../c/output/vars0_usa_a%d.csv'%a))
    for s in [0,1,2,3,4,5,6]:
        models_usa.append(pd.read_csv('../c/output/vars1_usa_t%d_s%d_c%d_r%d_d%d_a%d.csv'%(tau,s,c,r,d,a)))

    models_chn = []
    models_chn.append(pd.read_csv('../c/output/vars0_chn_a%d.csv'%a))
    for s in [0,1,2,3,4,5,6]:
        models_chn.append(pd.read_csv('../c/output/vars1_chn_t%d_s%d_c%d_r%d_d%d_a%d.csv'%(tau,s,c,r,d,a)))

    models_row = []
    models_row.append(pd.read_csv('../c/output/vars0_row_a%d.csv'%a))
    for s in [0,1,2,3,4,5,6]:
        models_row.append(pd.read_csv('../c/output/vars1_row_t%d_s%d_c%d_r%d_d%d_a%d.csv'%(tau,s,c,r,d,a)))
                        
    for m in models_usa:
        m['nx']=m['nx1']+m['nx2']
        m['lg']=m['l0']+m['l1']+m['l2']+m['l3']
        for s in [0,1,2,3,4]:
            m['nxs'+str(s)]=m['nxs'+str(s)+'-1']+m['nxs'+str(s)+'-2']
            m['exs'+str(s)]=m['exs'+str(s)+'-1']+m['exs'+str(s)+'-2']
            m['ims'+str(s)]=m['ims'+str(s)+'-1']+m['ims'+str(s)+'-2']
        m['nxsg'] = m['nxs0'] + m['nxs1'] + m['nxs2'] + m['nxs3']
        for col in m.columns:
            if('nx' in col):
                m[col] = m[col]/m.ngdp
    
    for m in models_chn:
        m['lg']=m['l0']+m['l1']+m['l2']+m['l3']
        m['nx']=m['nx0']+m['nx2']
        for s in [0,1,2,3,4]:
            m['nxs'+str(s)]=m['nxs'+str(s)+'-0']+m['nxs'+str(s)+'-2']
            m['exs'+str(s)]=m['exs'+str(s)+'-0']+m['exs'+str(s)+'-2']
            m['ims'+str(s)]=m['ims'+str(s)+'-0']+m['ims'+str(s)+'-2']
        m['nxsg'] = m['nxs0'] + m['nxs1'] + m['nxs2'] + m['nxs3']
        for col in m.columns:
            if('nx' in col):
                m[col] = m[col]/m.ngdp
    
    for m in models_row:
        m['lg']=m['l0']+m['l1']+m['l2']+m['l3']
        m['nx']=m['nx0']+m['nx1']
        for s in [0,1,2,3,4]:
            m['nxs'+str(s)]=m['nxs'+str(s)+'-0']+m['nxs'+str(s)+'-1']
            m['exs'+str(s)]=m['exs'+str(s)+'-0']+m['exs'+str(s)+'-1']
            m['ims'+str(s)]=m['ims'+str(s)+'-0']+m['ims'+str(s)+'-1']
        m['nxsg'] = m['nxs0'] + m['nxs1'] + m['nxs2'] + m['nxs3']
        for col in m.columns:
            if('nx' in col):
                m[col] = m[col]/m.ngdp
                

    return {'USA':models_usa,'CHN':models_chn,'ROW':models_row}

df = load_results(4)


##########################################################################
# make the plots

mperiods = df['USA'][0].period.values
slabels=['Up-hi','Up-lo','Dn-hi','Dn-lo','Total']

def pct_chg(df,country,variable,model_num):
    return (100*(df[country][model_num][variable]/df[country][0][variable] - 1)).values

def diff_pct_tot(df,country,variable,model_num):
    return (100*(df[country][model_num][variable]-df[country][0][variable])/df[country][0]['lg']).values

def pp_chg(df,country,variable,model_num):
    return (100*(df[country][model_num][variable]-df[country][0][variable])).values

# labor
for ts in [0,1,2,3,4,5,6]:
    
    fig,ax=plt.subplots(1,1,figsize=(5,5))
    cnt=0
    for s in [0,1,2,3]:
        data = diff_pct_tot(df,'USA','l%d'%s,ts+1)
        ax.plot(mperiods,data,color=colors[cnt],label=slabels[cnt])
        cnt=cnt+1

    data = diff_pct_tot(df,'USA','lg',ts+1)
    ax.plot(mperiods,data,color='black',label=slabels[-1],lw=2)    

    ax.legend(loc='best')
    fig.tight_layout()
    plt.savefig('output/labor_t%d.pdf'%(ts))
    plt.clf()

# prices
for ts in [0,1,2,3,4,5,6]:
    
    fig,ax=plt.subplots(1,1,figsize=(5,5))
    cnt=0
    for s in [0,1,2,3]:
        data = pct_chg(df,'USA','py%d'%s,ts+1)
        ax.plot(mperiods,data,color=colors[cnt],label=slabels[cnt])
        cnt=cnt+1

    ax.legend(loc='best')
    fig.tight_layout()
    plt.savefig('output/prices_t%d.pdf'%(ts))
    plt.clf()

# nx
for ts in [0,1,2,3,4,5,6]:
    
    fig,ax=plt.subplots(1,1,figsize=(5,5))
    cnt=0
    for s in [0,1,2,3]:
        data = pp_chg(df,'USA','nxs%d'%s,ts+1)
        ax.plot(mperiods,data,color=colors[cnt],label=slabels[cnt])
        cnt=cnt+1

    data = pp_chg(df,'USA','nxsg',ts+1)
    ax.plot(mperiods,data,color='black',label=slabels[-1],lw=2)    


    ax.legend(loc='best')
    fig.tight_layout()
    plt.savefig('output/nxs_t%d.pdf'%(ts))
    plt.clf()



