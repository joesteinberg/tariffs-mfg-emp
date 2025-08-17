# # make latex table
# trd=trd.groupby(['region','partner','sector']).mean().reset_index()

# sector_names = {'1-UPSTREAM':'Upstream goods',
#                 '2-DOWNSTREAM':'Downstream goods',
#                 '3-SERVICES':'Services',
#                 '4-CONSTRUCTION':'Construction',
#                 'TOT':'Total'}

# country_names = {'1-USA':'United States',
#                  '2-CHN':'China',
#                  '3-ROW':'Rest of world',
#                  'TOT':'Total'}
# partners = {'1-USA':['TOT','2-CHN','3-ROW'],'2-CHN':['TOT','1-USA','3-ROW'],'3-ROW':['TOT','1-USA','2-CHN']}
# panels = {'1-USA':'(a)','2-CHN':'(b)','3-ROW':'(c)'}

# def fmt_num(x):
#     if np.abs(x)<1.0e-6:
#         return '& -'
#     else:
#         return '& %0.2f'%x

# with open('output/icio_summary.tex','w') as file:
#     file.write('\\begin{tabular}{lccccc}\n')
#     file.write('\\toprule\n')
#     file.write('\\makecell{Quantity} & ')

#     file.write('\\makecell{Upstream\\\\goods} & ')
#     file.write('\\makecell{Downstream\\\\goods} &')
#     file.write('\\makecell{centering Services} &')
#     file.write('\\makecell{centering Construction} & ')
#     file.write('\\makecell{centering Total}\\\\\n')
#     file.write('\\midrule\n')

#     for c in ['1-USA','2-CHN','3-ROW']:
#         file.write('\\multicolumn{5}{l}{\\textit{'+panels[c]+' '+country_names[c]+'}}\\\\\n')
#         mask=trd.region==c

#         file.write('Value added')
#         mask2 = np.logical_and(mask,trd.partner=='TOT')
#         for s in sector_names.keys():
#             mask3=np.logical_and(mask2,trd.sector==s)
#             masked=trd[mask3]
#             val = 100.0*masked['VA']/masked['GDP']
#             file.write(fmt_num(val.iloc[0]))
#         file.write('\\\\\n')
        
#         for p in partners[c]:
#             mask2=np.logical_and(mask,trd.partner==p)
#             if p=='TOT':
#                 file.write('Exports')
#             else:
#                 file.write('\quad to ' + country_names[p])
#             for s in sector_names.keys():
#                 mask3=np.logical_and(mask2,trd.sector==s)
#                 masked=trd[mask3]
#                 val = 100.0*(masked['ex'])/masked['GDP']
#                 file.write(fmt_num(val.iloc[0]))
#             file.write('\\\\\n')

#         for p in partners[c]:
#             mask2=np.logical_and(mask,trd.partner==p)
#             if p=='TOT':
#                 file.write('Imports')
#             else:
#                 file.write('\quad from ' + country_names[p])
#             for s in sector_names.keys():
#                 mask3=np.logical_and(mask2,trd.sector==s)
#                 masked=trd[mask3]
#                 val = 100.0*(masked['im'])/masked['GDP']
#                 file.write(fmt_num(val.iloc[0]))
#             file.write('\\\\\n')

#         for p in partners[c]:
#             mask2=np.logical_and(mask,trd.partner==p)
#             if p=='TOT':
#                 file.write('Net exports')
#             else:
#                 file.write('\quad with ' + country_names[p])
#             for s in sector_names.keys():
#                 mask3=np.logical_and(mask2,trd.sector==s)
#                 masked=trd[mask3]
#                 val = 100.0*masked['tb']/masked['GDP']
#                 file.write(fmt_num(val.iloc[0]))
#             file.write('\\\\\n')

#         if(c!='ROW'):
#             file.write('\\\\\n')

#     file.write('\\bottomrule\n')
#     file.write('\\end{tabular}\n')

# figures
sector_names = {'1-UPSTREAM-HI':'Up-hi',
                '2-UPSTREAM-LO':'Up-lo',
                '3-DOWNSTREAM-HI':'Down-hi',
                '4-DOWNSTREAM-LO':'Down-lo',
                '5-SERVICES':'Svcs',
                '6-CONSTRUCTION':'Cnstr'}

country_names = {'1-USA':'United States',
                 '2-CHN':'China',
                 '3-ROW':'Rest of world',}

sectors = list(sector_names.keys())
snames = [sector_names[s] for s in sectors]
wrapped_snames= [ s.replace(' ', '\n') for s in snames ]

countries = list(country_names.keys())
cnames = [country_names[c] for c in countries]

partners = countries
pnames = cnames


df=trd.loc[trd.sector!='TOT',:].reset_index(drop=True)
df['trd']=df.ex+df.im
df['trd_M']=df.ex_M+df.im_M
df['trd2'] = df.trd
df['im2'] = df.im
df['ex2'] = df.ex
df['va2'] = df.VA

#print(df.ex_M)
for col in cols:
    if col in ['trd','trd_M','im_M','ex_M','im','ex']:
        df[col]=100*df[col]/df.VA
    else:
        df[col]=100*df[col]/df.GDP

#print(df.ex_M)

cols=['tb','ex','im','ex_M','im_M','trd2','trd']
avg_s = df[df.partner!='TOT'].groupby(['region','sector'])[cols].sum().reset_index()
#avg_s = df.groupby(['region','sector']).sum().reset_index()

avg_p = df[df.partner!='TOT'].groupby(['region','partner'])[cols].sum().reset_index()
#avg_va = df.groupby(['region','sector'])['va2'].mean().reset_index()

# figure 2: sectoral trade relative to VA
width=0.25
fig,axes=plt.subplots(2,2,figsize=(7,7),sharex='row',sharey=False)

# (a) sectoral gross imports/GDP
inds=np.arange(len(sectors[:-1]))

for i in range(len(countries)):
    data=np.zeros(len(inds))
    for j in range(len(sectors[:-1])):
        data[j]=avg_s.im[np.logical_and(avg_s.region==countries[i],
                                         avg_s.sector==sectors[j])].values[0]          
    if i==0:
        p1=axes[0,0].bar(inds,data,align='edge',width=width,color=colors[0],hatch=hatches[0],linewidth=1,alpha=0.99,edgecolor='black',)
    elif i==1:
        p2=axes[0,0].bar(inds+width,data,align='edge',width=width,color=colors[1],hatch=4*hatches[1],linewidth=1,alpha=0.99,edgecolor='black',)
    elif i==2:
        p3=axes[0,0].bar(inds+2*width,data,align='edge',width=width,color=colors[2],hatch=4*hatches[2],linewidth=1,alpha=0.99,edgecolor='black',)

axes[0,0].set_ylabel('percent sectoral value added')
axes[0,0].set_title('(a) Total imports',y=1.04,size=10)
axes[0,0].set_xticks(inds+width*1.5)
axes[0,0].set_xticklabels(wrapped_snames[:-1])
axes[0,0].set_xlim(-0.25,5.0)

# (b) sectoral intermediate imports/GDP
inds=np.arange(len(sectors[:-1]))

for i in range(len(countries)):
    data=np.zeros(len(inds))
    for j in range(len(sectors[:-1])):
        data[j]=avg_s.im_M[np.logical_and(avg_s.region==countries[i],
                                           avg_s.sector==sectors[j])].values[0]
          
    if i==0:
        p1=axes[0,1].bar(inds,data,align='edge',width=width,color=colors[0],hatch=hatches[0],linewidth=1,alpha=0.99,edgecolor='black',)
    elif i==1:
        p2=axes[0,1].bar(inds+width,data,align='edge',width=width,color=colors[1],hatch=4*hatches[1],linewidth=1,alpha=0.99,edgecolor='black',)
    elif i==2:
        p3=axes[0,1].bar(inds+2*width,data,align='edge',width=width,color=colors[2],hatch=4*hatches[2],linewidth=1,alpha=0.99,edgecolor='black',)

axes[0,1].set_title('(b) Intermediate imports',y=1.04,size=10)
axes[0,1].legend([p1,p2,p3],countries,loc='upper right',prop={'size':8})

# (c) sectoral gross imports/GDP
inds=np.arange(len(sectors[:-1]))

for i in range(len(countries)):
    data=np.zeros(len(inds))
    for j in range(len(sectors[:-1])):
        data[j]=avg_s.ex[np.logical_and(avg_s.region==countries[i],
                                         avg_s.sector==sectors[j])].values[0]          
    if i==0:
        p1=axes[1,0].bar(inds,data,align='edge',width=width,color=colors[0],hatch=hatches[0],linewidth=1,alpha=0.99,edgecolor='black',)
    elif i==1:
        p2=axes[1,0].bar(inds+width,data,align='edge',width=width,color=colors[1],hatch=4*hatches[1],linewidth=1,alpha=0.99,edgecolor='black',)
    elif i==2:
        p3=axes[1,0].bar(inds+2*width,data,align='edge',width=width,color=colors[2],hatch=4*hatches[2],linewidth=1,alpha=0.99,edgecolor='black',)

axes[1,0].set_ylabel('percent sectoral value added')
axes[1,0].set_title('(c) Total exports',y=1.04,size=10)
axes[1,0].set_xticks(inds+width*1.5)
axes[1,0].set_xticklabels(wrapped_snames[:-1],size=8)
axes[1,0].set_xlim(-0.25,5.0)

# (d) sectoral intermediate imports/GDP
inds=np.arange(len(sectors[:-1]))

for i in range(len(countries)):
    data=np.zeros(len(inds))
    for j in range(len(sectors[:-1])):
        data[j]=avg_s.ex_M[np.logical_and(avg_s.region==countries[i],
                                           avg_s.sector==sectors[j])].values[0]
          
    if i==0:
        p1=axes[1,1].bar(inds,data,align='edge',width=width,color=colors[0],hatch=hatches[0],linewidth=1,alpha=0.99,edgecolor='black',)
    elif i==1:
        p2=axes[1,1].bar(inds+width,data,align='edge',width=width,color=colors[1],hatch=4*hatches[1],linewidth=1,alpha=0.99,edgecolor='black',)
    elif i==2:
        p3=axes[1,1].bar(inds+2*width,data,align='edge',width=width,color=colors[2],hatch=4*hatches[2],linewidth=1,alpha=0.99,edgecolor='black',)

axes[1,1].set_title('(d) Intermediate exports',y=1.04,size=10)
axes[1,1].legend([p1,p2,p3],countries,loc='upper right',prop={'size':8})
fig.subplots_adjust(hspace=0.05,wspace=0.05)
fig.tight_layout()
plt.savefig('output/sectoral_trade.pdf')
plt.clf()

width=0.25


# # figure 3: net trade

vmax = lambda v: [max(x,0.0) for x in v]
vmin = lambda v: [min(x,0.0) for x in v]

def which_bottom(bottom_p, bottom_n, data):
    bottom=np.zeros(len(data))
    for i in range(len(data)):
        bottom[i]=bottom_p[i] if data[i]>0.0 else bottom_n[i]
    return bottom

fig,axes=plt.subplots(1,2,figsize=(7,3.5),sharex=True,sharey=True)

# (a) by partner
data=[np.zeros(3) for p in countries]
inds=range(len(countries))

for i in inds:
    p=partners[i]
    for j in inds[0:3]:
        c=countries[j]
        tmp=avg_p.tb[np.logical_and(avg_p.partner==p,avg_p.region==c)]
        if(len(tmp)>0):
            data[i][j]=tmp.values[0]

bottom_pos=np.zeros(3)
bottom_neg=np.zeros(3)
p1=axes[0].bar(inds[0:3],data[0],
               color=colors[0],hatch=hatches[0],linewidth=1,align='center',width=0.5,alpha=0.99,edgecolor='black',)
bottom_pos=bottom_pos+vmax(data[0])
bottom_neg=bottom_neg+vmin(data[0])

bottom=which_bottom(bottom_pos,bottom_neg,data[1])
p2=axes[0].bar(inds[0:3],data[1],bottom=bottom,
               color=colors[1],hatch=4*hatches[1],linewidth=1,align='center',width=0.5,alpha=0.99,edgecolor='black',)
bottom_pos=bottom_pos+vmax(data[1])
bottom_neg=bottom_neg+vmin(data[1])

bottom=which_bottom(bottom_pos,bottom_neg,data[2])
p3=axes[0].bar(inds[0:3],data[2],bottom=bottom,
               color=colors[2],hatch=4*hatches[2],linewidth=1,align='center',width=0.5,alpha=0.99,edgecolor='black',)
bottom_pos=bottom_pos+vmax(data[2])
bottom_neg=bottom_neg+vmin(data[2])


axes[0].plot(range(-1,4),np.zeros(len(range(-1,4))),color='black',linestyle='-')
axes[0].set_xticks(inds)
axes[0].set_xticklabels(countries)
#axes[0].set_yticks([-3,0,3,6,9,12])
#axes[0].set_ylim(-3,12)
axes[0].set_xlim(-0.5,2.5)
axes[0].set_ylabel('percent GDP')
axes[0].set_title('(a) By partner',y=1.04,size=10)
axes[0].legend([p1,p2,p3],cnames,loc='upper left',prop={'size':8},ncol=1)

# (b) by sector
data=[np.zeros(3) for s in sectors]
inds=range(len(countries))

for i in range(len(sectors)):
    s=sectors[i]
    for j in inds:
        c=countries[j]
        tmp=avg_s.tb[np.logical_and(avg_s.sector==s,avg_s.region==c)]
        if(len(tmp)>0):
            data[i][j]=tmp.values[0]

bottom_pos=np.zeros(3)
bottom_neg=np.zeros(3)
p1=axes[1].bar(inds,data[0],color=colors[0],hatch=hatches[0],align='center',width=0.5,alpha=0.99,edgecolor='black',)
bottom_pos=bottom_pos+vmax(data[0])
bottom_neg=bottom_neg+vmin(data[0])

bottom=which_bottom(bottom_pos,bottom_neg,data[1])
p2=axes[1].bar(inds,data[1],bottom=bottom,
               color=colors[1],hatch=4*hatches[1],linewidth=1,align='center',width=0.5,alpha=0.99,edgecolor='black',)
bottom_pos=bottom_pos+vmax(data[1])
bottom_neg=bottom_neg+vmin(data[1])

bottom=which_bottom(bottom_pos,bottom_neg,data[2])
p3=axes[1].bar(inds,data[2],bottom=bottom,
               color=colors[2],hatch=4*hatches[2],linewidth=1,align='center',width=0.5,alpha=0.99,edgecolor='black',)
bottom_pos=bottom_pos+vmax(data[2])
bottom_neg=bottom_neg+vmin(data[2])

bottom=which_bottom(bottom_pos,bottom_neg,data[3])
p4=axes[1].bar(inds,data[3],bottom=bottom,
               color=colors[3],hatch=4*hatches[3],linewidth=1,align='center',width=0.5,alpha=0.99,edgecolor='black',)
bottom_pos=bottom_pos+vmax(data[3])
bottom_neg=bottom_neg+vmin(data[3])

bottom=which_bottom(bottom_pos,bottom_neg,data[4])
p5=axes[1].bar(inds,data[4],bottom=bottom,
               color=colors[4],hatch=4*hatches[4],linewidth=1,align='center',width=0.5,alpha=0.99,edgecolor='black',)
bottom_pos=bottom_pos+vmax(data[4])
bottom_neg=bottom_neg+vmin(data[4])


axes[1].plot(range(-1,4),np.zeros(len(range(-1,4))),color='black',linestyle='-')
axes[1].set_title('(b) By sector',y=1.04,size=10)
axes[1].legend([p1,p2,p3,p4,p5],snames,loc='upper left',prop={'size':8},ncol=2)

fig.subplots_adjust(hspace=0.05,wspace=0.05)
fig.tight_layout()
plt.savefig('output/trade_balances.pdf')
plt.clf()





plt.close('all')
