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

colors=['#e41a1c','#377eb8','#984ea3','#4daf4a','#ff7f00','#ffff33']
hatches=[None,'x','/','.','+']

vmax = lambda v: [max(x,0.0) for x in v]
vmin = lambda v: [min(x,0.0) for x in v]

def which_bottom(bottom_p, bottom_n, data):
    bottom=np.zeros(len(data))
    for i in range(len(data)):
        bottom[i]=bottom_p[i] if data[i]>0.0 else bottom_n[i]
    return bottom

from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler

#############################################################################
print('Processing the raw data...')

# load the CSV in matrix format
csv = pd.read_csv('../../data/2020_SML.csv').\
    rename(columns={'V1':'row_label'})

# melt it to long format
melted = pd.melt(csv,id_vars = 'row_label',var_name='col_label',
                 value_vars = csv.columns[1:])


# assign country and industry codes
melted['row_country'] = ''
melted['col_country'] = ''
melted['row_ind'] = ''
melted['col_ind'] = ''

mask = melted.row_label.isin(['TLS','VA','OUT'])
melted.loc[mask,'row_country'] = 'TOT'
melted.loc[mask,'row_ind'] = melted.loc[mask,'row_label']
melted.loc[~mask,'row_country'] = melted.loc[~mask,'row_label'].str[0:3]
melted.loc[~mask,'row_ind'] = melted.loc[~mask,'row_label'].str[4:]

mask = melted.col_label.isin(['OUT'])
melted.loc[mask,'col_country'] = 'TOT'
melted.loc[mask,'col_ind'] = melted.loc[mask,'col_label']
melted.loc[~mask,'col_country'] = melted.loc[~mask,'col_label'].str[0:3]
melted.loc[~mask,'col_ind'] = melted.loc[~mask,'col_label'].str[4:]

# assign use types to columns
melted['col_use'] = ''

def which_use(ind):
    if ind in(['HFCE','NPISH','GGFC','GFCF','INVNT','DPABR']):
        return 'FIN'
    elif ind=='OUT':
        return 'TOT'
    else:
        return 'INT'
    
melted.loc[:,'col_use'] = melted.loc[:,'col_ind'].apply(which_use)


# assign region aggregation
regions = ['USA','CHN']
def which_region(country):
    if country=='TOT' or country in regions:
        return country
    else:
        return 'ROW'

melted['row_region'] = melted.row_country.apply(which_region)
melted['col_region'] = melted.col_country.apply(which_region)

melted.loc[ (melted.col_ind=='INVNT') & (melted.value<0),'value'] = 0


# aggregate non-mfg sectors
def which_sector(ind):
    if ind in ['TOT','OUT']:
        return 'GO'
    elif ind in ['VA','TLS']:
        return 'VA'
    elif ind in ['HFCE','NPISH','GGFC','DPABR']:
        return 'CC'
    elif ind in ['GFCF','INVNT']:
        return 'INV'
    elif ind[0] in ['A','B','C']:
        return ind
    elif ind=='F':
        return 'CONS'
    else:
        return 'SVCS'
    
melted['row_sector'] = melted.row_ind.apply(which_sector)
melted['col_sector'] = melted.col_ind.apply(which_sector)


#############################################################################
print('\nComputing upstreamness in goods industries...')

# construct IO matrix aggregated across regions
agged = melted.groupby(['row_region','row_sector','col_use','col_region','col_sector'])['value'].sum().reset_index()
agged = agged.sort_values(by=['row_region','row_sector','col_use','col_region','col_sector'])

# compute direct requirement coefficients
wide = agged.loc[agged.row_region!='TOT',:].pivot_table(values='value',index=['row_region','row_sector'],columns=['col_use','col_region','col_sector'])

Y = wide['TOT']
M = wide['INT']
D = M.copy(deep=True)
D.values[:,:] = M.values / Y.values.transpose()

# compute upststreamness measure from Antras et al. 2012
# ''Measuring the Upstreamness of Production and Trade Flows''
# https://scholar.harvard.edu/files/antras/files/acfh_published.pdf
# U_i... equals the dollar amount by which 
# output of all sectors increases following a one 
# dollar increase in value added in sector i. This is 
# a standard measure of cost-push effects or total 
# forward linkages in supply-side I-O models and 
# is intuitively increasing in upstreamness


Delta = D.copy(deep=True)
Delta.values[:,:] = D.values[:,:] * Y.values.transpose() / Y.values
U = Y.copy(deep=True)
U.values[:] = np.matmul(np.linalg.inv(np.eye(U.shape[0]) - Delta), np.ones((U.shape[0],1)))
U.columns = U.columns.droplevel(0)
U.columns = ['U']
U = U.reset_index()
U = U.loc[~(U.row_sector.isin(['AGRI','MINE','SVCS','CONS'])),:]

# display results
U = U.rename(columns={'row_sector':'ind'})
names = pd.read_csv('../../data/industry_names_elasts.csv')[['ind','name']]
U = pd.merge(left=U,right=names,how='left',on='ind')
U2 = U.groupby(['ind','name'])['U'].mean().sort_values(ascending=False).reset_index()
#U2['upstream'] = np.where(U2.U > U2.U.median(),True,False)
print(U2[['ind','name','U']].sort_values(by='U'))


#############################################################################
print('\nComputing trade elasticities...')
elasts = pd.read_csv('../../data/industry_names_elasts.csv')[['ind','elast_CP']].rename(columns={'ind':'row_sector'})

merged = pd.merge(left=agged,right=elasts,how='left',on='row_sector')
trd = merged[ (merged.col_region!=merged.row_region) &
              (merged.row_region!='TOT') &
              (merged.col_region!='TOT') &
              (merged.col_use!='TOT') &
              (merged.row_sector!='TOT') &
              (merged.row_sector!='VA') ]

# first, compute weighted average elasticity by industry
wavg = lambda x: np.average(x,weights=trd.loc[x.index,'value'])
e = trd.groupby(['row_sector'])[['elast_CP']].agg(wavg).reset_index()

#############################################################################
print('\nDefining sectoral aggregation using hierarchical clustering...')

scheme = e.loc[:,['row_sector','elast_CP']]
scheme = pd.merge(left=scheme,right=U2.rename(columns={'ind':'row_sector'}),how='left',on='row_sector')
scheme2 = scheme.loc[~(scheme.name.isnull()),['row_sector','name','elast_CP','U']]

tmp = np.log(scheme2.elast_CP)
X = np.zeros((len(scheme2),2))
X[:,0] = scheme2.U
X[:,1] = tmp

# Cluster on the first dimension
X_dim1 = X[:,0].reshape(-1, 1)
    
# Scale the data if necessary 
scaler_dim1 = StandardScaler()
X_dim1_scaled = scaler_dim1.fit_transform(X_dim1)

# Apply Agglomerative Clustering
n_clusters_dim1 = 2
clustering_dim1 = AgglomerativeClustering(n_clusters=n_clusters_dim1)
labels_dim1 = clustering_dim1.fit_predict(X_dim1_scaled)

# Cluster on the second dimension within each cluster from the first dimension
final_labels = np.zeros(len(X), dtype=int)
current_cluster_id = 0

for i in range(n_clusters_dim1):
    # Get data points belonging to the current cluster from the first dimension
    cluster_indices = np.where(labels_dim1 == i)[0]
    X_dim2_subset = X[cluster_indices,1].reshape(-1, 1)

    if len(X_dim2_subset) > 0:
        scaler_dim2 = StandardScaler()
        X_dim2_subset_scaled = scaler_dim2.fit_transform(X_dim2_subset)
        
        # Apply Agglomerative Clustering
        n_clusters_dim2 = 2
        clustering_dim2 = AgglomerativeClustering(n_clusters=n_clusters_dim2)
        labels_dim2_subset = clustering_dim2.fit_predict(X_dim2_subset_scaled)

        # Assign new, unique cluster IDs
        for j, original_index in enumerate(cluster_indices):
            final_labels[original_index] = current_cluster_id + labels_dim2_subset[j]
        current_cluster_id += n_clusters_dim2

scheme2['row_sector_num'] = final_labels

def which_sector_name(n):
    if n==0:
        return 'UP-HI'
    elif n==1:
        return 'UP-LO'
    elif n==2:
        return 'DN-LO'
    elif n==3:
        return 'DN-HI'

scheme2['row_sector3'] = scheme2.row_sector_num.apply(which_sector_name)

# print the resulting classification
print(scheme2.sort_values(by='row_sector3'))

# compute weighted average elasticities for each group
merged = pd.merge(left=agged,right=scheme2[['row_sector','row_sector3','elast_CP']],how='left',on='row_sector')
                   
trd = merged[ (merged.col_region!=merged.row_region) &
              (merged.row_region!='TOT') &
              (merged.col_region!='TOT') &
              (merged.col_use!='TOT') &
              (merged.row_sector!='TOT') &
              (merged.row_sector!='VA') ]

wavg = lambda x: np.average(x,weights=trd.loc[x.index,'value'])
e2 = trd.groupby(['row_sector3'])[['elast_CP']].agg(wavg).reset_index()
print(e2.loc[~(e2.elast_CP.isnull()),:])

# merge the groups back on to the main dataset
# UP-HI = OIL
# UP-LO = STEEL
# DN-HI = CONSUMER GOODS
# DN-LO = MACHINERY & AUTOS
agged = pd.merge(left=agged,right=scheme2[['row_sector','row_sector3']],how='left',on='row_sector')
agged = pd.merge(left=agged,right=scheme2[['row_sector','row_sector3']].rename(columns={'row_sector':'col_sector','row_sector3':'col_sector3'}),how='left',on='col_sector')
agged.loc[agged.row_sector3.isnull(),'row_sector3'] = agged.loc[agged.row_sector3.isnull(),'row_sector']
agged.loc[agged.col_sector3.isnull(),'col_sector3'] = agged.loc[agged.col_sector3.isnull(),'col_sector']

#############################################################################
print('\nComputing sectoral employment shares...')

lc = pd.read_excel('../../data/ComponentsOfVa.xlsx',sheet_name='UVCT2-A',skiprows=7).rename(columns={'ICIO Industry Code':'ind','2023':'lcomp'})
lc = lc.loc[~(lc.ind.isnull()),['ind','lcomp']].groupby('ind').sum().reset_index()
lc['lshare'] = lc.lcomp/lc.lcomp.sum()

icio_inds = list(scheme.row_sector)
def find_sector(bea_ind):
    if bea_ind == 'F':
        return 'CONS'
    else:
        for i in icio_inds:
            if i in bea_ind:
                return scheme2.loc[scheme.row_sector==i,'row_sector3'].values[0]
        return 'SVCS'

lc['sector'] = lc.ind.apply(find_sector)
lc2=lc.groupby('sector').lshare.sum()
ltot = lc2['UP-HI'] + lc2['UP-LO'] + lc2['DN-HI'] + lc2['DN-LO']

print(lc2)

#############################################################################
print('\nCreating sector summary figure + table...')

# write table with names
with open('output/sectors.tex','w') as file:
        
    file.write('\\begin{tabular}{lp{6cm}cccc}')
    file.write('\\toprule\n')
    file.write('Sector & Industries & Upstreamness & \\makecell{Trade\\\\elasticity} & \\makecell{Share of\\\\mfg. emp.}\\\\\n')
    file.write('\\midrule\n')

    file.write("\\textcolor[HTML]{%s}{``Oil''} & "%(colors[0][1:]))
    cnt=0
    names=scheme2.loc[scheme2.row_sector_num==0,'name'].values
    for n in names:
        cnt += 1
        file.write('%s'%n)
        if(cnt<len(names)):
           file.write(', ')
           
    file.write('& %0.1f' % (scheme2.loc[scheme2.row_sector_num==0,'U'].mean()))
    file.write('& %0.1f' % (e2.loc[e2.row_sector3=='UP-HI','elast_CP'].values[0]))
    file.write('& %0.1f\\\\\n' % (100*lc2['UP-HI']/ltot))

    file.write("\\textcolor[HTML]{%s}{``Steel''} & "%(colors[1][1:]))
    cnt=0
    names=scheme2.loc[scheme2.row_sector_num==1,'name'].values
    for n in names:
        cnt += 1
        file.write('%s'%n)
        if(cnt<len(names)):
           file.write(', ')

    file.write('& %0.1f' % (scheme2.loc[scheme2.row_sector_num==1,'U'].mean()))
    file.write('& %0.1f' % (e2.loc[e2.row_sector3=='UP-LO','elast_CP'].values[0]))
    file.write('& %0.1f\\\\\n' % (100*lc2['UP-LO']/ltot))

    file.write("\\textcolor[HTML]{%s}{``Toys''} & "%(colors[3][1:]))
    cnt=0
    names=scheme2.loc[scheme2.row_sector_num==3,'name'].values
    for n in names:
        cnt += 1
        file.write('%s'%n)
        if(cnt<len(names)):
           file.write(', ')

    file.write('& %0.1f' % (scheme2.loc[scheme2.row_sector_num==3,'U'].mean()))
    file.write('& %0.1f' % (e2.loc[e2.row_sector3=='DN-HI','elast_CP'].values[0]))
    file.write('& %0.1f\\\\\n' % (100*lc2['DN-HI']/ltot))

    file.write("\\textcolor[HTML]{%s}{``Cars''} & "%(colors[2][1:]))
    cnt=0
    names=scheme2.loc[scheme2.row_sector_num==2,'name'].values
    for n in names:
        cnt += 1
        file.write('%s'%n)
        if(cnt<len(names)):
           file.write(', ')

    file.write('& %0.1f' % (scheme2.loc[scheme2.row_sector_num==2,'U'].mean()))
    file.write('& %0.1f' % (e2.loc[e2.row_sector3=='DN-LO','elast_CP'].values[0]))
    file.write('& %0.1f\\\\\n' % (100*lc2['DN-LO']/ltot))

    file.write('\\bottomrule\n')
    file.write('\\end{tabular}')

# figure
scheme2['log_e'] = np.log(scheme2.elast_CP)
scheme3 = pd.merge(left=scheme2,right=lc.rename(columns={'ind':'row_sector'})[['row_sector','lshare']],how='left',on='row_sector')

final_labels2 = final_labels.copy()
for i in range(len(final_labels)):
    if final_labels[i]==2:
        final_labels2[i]=3
    elif final_labels[i]==3:
        final_labels2[i]=2

# figure
fig, ax = plt.subplots(figsize=(3.25,3.4),tight_layout = {'pad': 0})
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)
ax.yaxis.label.set_size(12)
sns.despine()
ax.scatter(scheme3.U,scheme3.log_e, c=[colors[k] for k in final_labels2])

tmpv = 0.5*scheme3[scheme3.row_sector_num<2]['U'].min() +  0.5*scheme3[scheme3.row_sector_num>=2]['U'].max()
ax.axvline(tmpv,color='black',lw=1)

tmph = np.log(0.5*scheme3.loc[scheme3.row_sector_num==3,'elast_CP'].min() + 0.5*scheme3.loc[scheme3.row_sector_num==2,'elast_CP'].max())
ax.plot([scheme3.U.min(),tmpv],[tmph,tmph],color='black',lw=1)

tmph = np.log(0.5*scheme3.loc[scheme3.row_sector_num==0,'elast_CP'].min() + 0.5*scheme3.loc[scheme3.row_sector_num==1,'elast_CP'].max())
ax.plot([tmpv,scheme3.U.max()],[tmph,tmph],color='black',lw=1)

yl = ax.get_yticks()
ax.set_yticklabels([ '%0.0f'%(np.exp(y)) for y in yl])
ax.set_xlabel('Upstreamness')
ax.set_ylabel('Trade elasticity')
fig.tight_layout()
plt.savefig('output/fig_data_upstream_vs_elast.pdf')

#############################################################################
print('\nAggregating IO matrix across regions and sectors...')

# now we can overwrite the aggregated sector
agged['row_sector'] = agged['row_sector3']
agged['col_sector'] = agged['col_sector3']

agged2 = agged.groupby(['row_region','row_sector','col_use','col_region','col_sector'])['value'].sum().reset_index()
agged2 = agged2.sort_values(by=['row_region','row_sector','col_use','col_region','col_sector'])

# order region and sectors as desired
agged2.loc[agged2.row_region=='USA','row_region'] = '1-USA'
agged2.loc[agged2.row_region=='CHN','row_region'] = '2-CHN'
agged2.loc[agged2.row_region=='ROW','row_region'] = '3-ROW'
agged2.loc[agged2.col_region=='USA','col_region'] = '1-USA'
agged2.loc[agged2.col_region=='CHN','col_region'] = '2-CHN'
agged2.loc[agged2.col_region=='ROW','col_region'] = '3-ROW'

agged2.loc[agged2.row_sector=='UP-HI','row_sector'] = '1-UPSTREAM-HI'
agged2.loc[agged2.row_sector=='UP-LO','row_sector'] = '2-UPSTREAM-LO'
agged2.loc[agged2.row_sector=='DN-HI','row_sector'] = '3-DOWNSTREAM-HI'
agged2.loc[agged2.row_sector=='DN-LO','row_sector'] = '4-DOWNSTREAM-LO'
agged2.loc[agged2.row_sector=='SVCS','row_sector'] = '5-SERVICES'
agged2.loc[agged2.row_sector=='CONS','row_sector'] = '6-CONSTRUCTION'

agged2.loc[agged2.col_sector=='UP-HI','col_sector'] = '1-UPSTREAM-HI'
agged2.loc[agged2.col_sector=='UP-LO','col_sector'] = '2-UPSTREAM-LO'
agged2.loc[agged2.col_sector=='DN-HI','col_sector'] = '3-DOWNSTREAM-HI'
agged2.loc[agged2.col_sector=='DN-LO','col_sector'] = '4-DOWNSTREAM-LO'
agged2.loc[agged2.col_sector=='SVCS','col_sector'] = '5-SERVICES'
agged2.loc[agged2.col_sector=='CONS','col_sector'] = '6-CONSTRUCTION'

# separate into main components in same structure as in NAFTA paper
intermediates = agged2[(agged2.col_region!='TOT') &
                       (agged2.row_region!='TOT') &
                       (agged2.col_use=='INT')].groupby(['col_region','col_sector','row_region','row_sector'])['value']\
                       .sum().reset_index().rename(columns={'value':'M'})

consumption = agged2[(agged2.col_sector=='CC') &
                     (agged2.col_region!='TOT') &
                     (agged2.row_region!='TOT')].groupby(['col_region','row_region','row_sector'])['value']\
                     .sum().reset_index().rename(columns={'value':'C'})

investment = agged2[(agged2.col_sector=='INV') &
                    (agged2.col_region!='TOT') &
                    (agged2.row_region!='TOT')].groupby(['col_region','row_region','row_sector'])['value']\
                    .sum().reset_index().rename(columns={'value':'I'})

value_added = agged2[(agged2.col_region!='TOT') &
                     (agged2.col_sector!='TOT') &
                     (agged2.col_use=='INT') &
                     (agged2.row_sector=='VA')].groupby(['col_region','col_sector'])['value']\
                     .sum().reset_index().rename(columns={'value':'VA'})

gross_output = agged2[(agged2.col_region!='TOT') &
                      (agged2.col_sector!='TOT') &
                      (agged2.col_use=='INT') &
                      (agged2.row_sector=='GO')].groupby(['col_region','col_sector'])['value']\
                      .sum().reset_index().rename(columns={'value':'GO'})

final_demand = pd.merge(left=consumption,right=investment,
                        how='left',
                        on=['col_region','row_region','row_sector'])

output = pd.merge(left=value_added,right=gross_output,
                  how='left',
                  on=['col_region','col_sector'])

#############################################################################
print('\nChecking market clearing...')


# check market clearing
msums = intermediates.groupby(['row_region','row_sector'])['M'].sum().reset_index()
msums.rename(columns={'row_region':'region','row_sector':'sector'},inplace=True)

fsums = final_demand.groupby(['row_region','row_sector'])[['C','I']].sum().reset_index()
fsums.rename(columns={'row_region':'region','row_sector':'sector'},inplace=True)

gsums = output[['col_region','col_sector','GO']]
gsums = gsums.rename(columns={'col_region':'region','col_sector':'sector'})

sums = pd.merge(left=msums,right=fsums,how='left',on=['region','sector'])
sums = pd.merge(left=sums,right=gsums,how='left',on=['region','sector'])
sums['diff'] = (sums.GO - sums.M - sums.C - sums.I)
sums['diff'] = sums['diff']/sums['GO']
test = sum(sums['diff']>1e-4)

if test>0:
    print('Market clearing failure!')

#############################################################################
print('\nAppyling assumptions and ensuring IO matrix is balanced...')


# apply assumption: construction is purely nontraded and used only for consumption
final_demand.loc[ (final_demand.row_sector=='6-CONSTRUCTION') & (final_demand.col_region != final_demand.row_region), 'I']=0
intermediates.loc[ (intermediates.row_sector=='6-CONSTRUCTION') , 'M']=0
final_demand.loc[ (final_demand.row_sector=='6-CONSTRUCTION'), 'C']=0
    
# ensure IO matrix is balanced
nc  = len(final_demand.row_region.unique())
ns = len(final_demand.row_sector.unique())

rowsums = np.zeros( nc*ns + 1 )
colsums = np.zeros( nc*ns + nc*2 )

MM = intermediates.pivot_table(values='M', index=['row_region','row_sector'], columns=['col_region','col_sector'])
VV = output['VA'].values.reshape((1,nc*ns))
FF = final_demand.pivot_table(values=['C','I'], index=['row_region','row_sector'], columns=['col_region'])
VV = np.hstack((VV,np.zeros((1,nc*2))))

iomat=np.vstack( ( np.hstack((MM,FF)) , VV ) )

for row in range(0,nc*ns + 1):
    rowsums[row] = np.sum(iomat[row,:])

for col in range(0,nc*ns + nc*2):
    colsums[col] = np.sum(iomat[:,col])

def coeffs(iomat):
    # Given world IO matrix (iomat), calculates IO coefficients and returs them in A

    A=np.zeros(iomat.shape)
    for col in range(0,A.shape[1]):
        A[:,col] = iomat[:,col]/np.sum(iomat[:,col])
    return A

def ras(iomat0,rowsums1,colsums1):
    # Given an initial IO matrix (iomat), and desired rowsums (rowsums1) and colsums (colsums1),
    # performs the RAS balancing procedure. Returns a new IO matrix (iomat) that is consistent
    # with desired row- and colsums.

    A0 = coeffs(iomat0)
    iomat = np.dot(A0,np.diag(colsums1))

    go=True
    iter=0
    maxit=10000
    tol=1.0e-8

    while go:
        iter=iter+1
        rowsums = np.sum(iomat,axis=1)
        r = np.divide(rowsums1,rowsums)
        iomat = np.dot(np.diag(r),iomat)
        colsums = np.sum(iomat,axis=0)
        s = np.divide(colsums1,colsums)
        iomat = np.dot(iomat,np.diag(s))
        colsums = np.sum(iomat,axis=0)
        rowsums = np.sum(iomat,axis=1)

        norm1 = max(np.divide(abs(rowsums-rowsums1),rowsums1))
        norm2 = max(np.divide(abs(colsums-colsums1),colsums1))
        if((norm1 <tol and norm2 <tol) or iter == maxit):
            go=False

    if iter==maxit:
        print('\tRAS iteration did not converge!')
        print('\titer = ', iter, ' diff = ', max(norm1,norm2))
    else:
        print('\tRAS converged after ',str(iter),' iterations')


    return iomat

colsums[0:(nc*ns)] = rowsums[0:(nc*ns)] # make sure markets clear: gross output = total demand for each country/sector
rowsums[-1] = colsums[(nc*ns):].sum() # world value added must equal world final demand
iomat2 = ras(iomat,rowsums,colsums) # run RAS


#############################################################################
print('Writing output...')

countries = ['USA','CHN','ROW']
sectors = ['$G_{UH}$','$G_{UL}$','$G_{DH}$','$G_{DL}$','$S$','$H$']

def write_iomat_csv(iomat,fname):
    usgdp = iomat[-1,0:ns].sum()
    iomat2 = np.vstack((iomat,np.sum(iomat,axis=0).reshape((1,nc*ns+nc*2))))
    iomat2 = np.hstack((iomat2,np.sum(iomat2,axis=1).reshape((nc*ns+2,1))))
    iomat2 = 100*iomat2/usgdp
    np.savetxt(fname=fname+'.txt',X=iomat2,fmt='%0.15f',delimiter=' ')

def write_iomat_latex(iomat,rowsums,colsums,fname):
    usgdp = iomat[-1,0:ns].sum()
    iomat2 = 100*iomat[:,:]/usgdp
    rowsums2 = 100*rowsums/usgdp
    colsums2 = 100*colsums/usgdp

    M=iomat2[0:(nc*ns),0:(nc*ns)]
    V=iomat2[-1,0:(nc*ns)]
    Fc=iomat2[0:(nc*ns),(nc*ns):+((nc*ns)+nc)]
    Fx=iomat2[0:(nc*ns),((nc*ns)+nc):]

    with open(fname + '.tex','w') as file:
        
        file.write('\\begin{tabular}{cc ')
        for i in range(0,nc*ns):
            file.write('c')
        file.write(' ')
        for i in range(0,nc*2):
            file.write('c')
        file.write('}\n')
        file.write('\\toprule\n')

        file.write('&&\\multicolumn{'+str(ns*nc)+'}{c}{Intermediate inputs}')
        file.write('&\\multicolumn{'+str(nc*2)+'}{c}{Final demand}\\\\\n')

        x=3
        file.write('\\cmidrule(rl){'+str(x)+'-'+str(x+ns*nc-1)+'}')
        x = x+ns*nc
        file.write('\\cmidrule(rl){'+str(x)+'-'+str(x+nc*2-1)+'}\n')
        
        # country names (intermediate input portion)
        file.write('&')
        for c in countries:
            file.write('& \\multicolumn{'+str(ns)+'}{c}{'+c+'}')
        #file.write('&')
        
        # country names (final demand portion)
        for c in countries:
            file.write('& \\multicolumn{2}{c}{'+c+'}')
        file.write('\\\\\n')

        # underline country names (intermediate input portion)
        x=3
        for i in range(nc):
            file.write('\\cmidrule(rl){'+str(x)+'-'+str(x+ns-1)+'}')
            x=x+ns

        # underline country names (final demand portion)
        for i in range(nc):
            file.write('\\cmidrule(rl){'+str(x)+'-'+str(x+1)+'}')
            x=x+2
        file.write('\n')

        # sector names  (intermediate input portion)
        file.write('&')
        for c in countries:
            for s in sectors:
                file.write('&' + s)

        # sector names (final demand portion)
        for c in countries:
            file.write('& $C$ & $X$')
        file.write('\\\\\n')
        file.write('\\midrule\n')

        # main data rows
        for i in range(0,nc):
            
            file.write('\\multirow{'+str(ns)+'}{*}{\\begin{sideways}'+countries[i]+'\\end{sideways}}')
            
            for ii in range(0,ns):
                
                file.write('&'+sectors[ii])

                # intermediate inputs
                for j in range(0,nc):                    
                    for jj in range(0,ns):
                        tmpstr = '-'
                        if M[i*ns+ii][j*ns+jj] > 1e-6:
                            tmpstr = locale.format_string('%0.2f',M[i*ns+ii,j*ns+jj],grouping=True)
                        file.write('&'+tmpstr)

                # final demand
                for j in range(0,nc):
                    tmpstr='-'
                    if Fc[i*ns+ii,j]>1e-6:
                        tmpstr = locale.format_string('%0.2f',Fc[i*ns+ii,j],grouping=True)
                    file.write('&'+tmpstr)

                    tmpstr='-'
                    if abs(Fx[i*ns+ii,j]>1e-6):
                        tmpstr = locale.format_string('%0.2f',Fx[i*ns+ii,j],grouping=True)
                    file.write('&'+tmpstr)
                

                file.write('\\\\\n')
            file.write('\\midrule\n')
            

        # value added row
        file.write('& VA')
        for i in range(0,nc):
            for ii in range(0,ns):
                tmpstr='-'
                if V[i*ns+ii]>1e-6:
                    tmpstr = locale.format_string('%0.2f',V[i*ns+ii],grouping=True)
                file.write('&'+tmpstr)
        for x in range(0,nc*2):
            file.write('&--')

        file.write('\\\\\n')

        # gross output row
        file.write('\\midrule\n')
        file.write('& GO')
        for i in range(0,nc):
            for ii in range(0,ns):
                tmpstr='-'
                if colsums2[i*ns+ii]>1e-6:
                    tmpstr = locale.format_string('%0.2f',colsums2[i*ns+ii],grouping=True)
                file.write('&'+tmpstr)
        for x in range(0,nc*2):
            file.write('&--')

        file.write('\\\\\n')
            
        file.write('\\bottomrule\n')
        file.write('\\end{tabular}\n')
        #file.write('\\end{center}\n')
        #file.write('\\end{table}\n')
        #file.write('\\end{landscape}\n')


write_iomat_csv(iomat2,'output/iomat')
write_iomat_latex(iomat2,rowsums,colsums,'output/iomat')

##################################################################################
# descriptive tables/figures
print('Figures summarizing trade flows')

def create_fig():
    fig, ax = plt.subplots(figsize=(3.25,3.25),tight_layout = {'pad': 0})
    ax.tick_params(axis='both', labelsize=10)
    ax.yaxis.label.set_size(10)
    sns.despine()
    return fig,ax

# intermediate trade
m_trd =  intermediates.groupby(['col_region','row_region','row_sector'])['M'].sum().reset_index()
m_trd=m_trd[m_trd.col_region != m_trd.row_region]
im_m = m_trd.rename(columns={'col_region':'region','row_sector':'sector','row_region':'partner','M':'im_M'})
ex_m = m_trd.rename(columns={'row_region':'region','row_sector':'sector','col_region':'partner','M':'ex_M'})
m_trd2 = pd.merge(left=ex_m,right=im_m,how='left',on=['region','partner','sector'])

# final trade
f_trd = final_demand[final_demand.col_region != final_demand.row_region]
im_f = f_trd.rename(columns={'col_region':'region','row_sector':'sector','row_region':'partner','C':'im_C','I':'im_I'})
ex_f = f_trd.rename(columns={'row_region':'region','row_sector':'sector','col_region':'partner','C':'ex_C','I':'ex_I'})
f_trd2 = pd.merge(left=ex_f,right=im_f,how='left',on=['region','partner','sector'])

# merge and calculate totals + balances
trd = pd.merge(left=m_trd2,right=f_trd2,how='left',on=['region','partner','sector'])

for d in ['ex','im']:
    trd[d+'_F'] = trd[d+'_C']+trd[d+'_I']
    trd[d] = trd[d+'_M']+trd[d+'_F']

for u  in ['_M','_C','_I','_F','']:
    trd['tb'+u] = trd['ex'+u] - trd['im'+u]

# aggregate by sector and append
cols = []
for d in ['ex','im','tb']:
    for u in ['','_M','_F','_C','_I']:
        cols.append(d+u)

g = trd.groupby(['region','partner'])
sums = g[cols].sum().reset_index()
sums['sector'] = 'TOT'
trd = pd.concat([trd,sums])

# aggregate by country and append
g = trd.groupby(['region','sector'])
sums = g[cols].sum().reset_index()
sums['partner'] = 'TOT'
trd = pd.concat([trd,sums])
trd = trd.sort_values(['region','partner','sector']).reset_index(drop=True)

# merge on value added
go = output.groupby(['col_region','col_sector'])['GO'].sum().reset_index()
go.rename(columns={'col_region':'region','col_sector':'sector'},inplace=True)
trd = pd.merge(left=trd,right=go,how='left',on=['region','sector'])

# merge on consumption
cons = final_demand.groupby(['col_region','row_sector'])['C'].sum().reset_index()
cons.rename(columns={'col_region':'region','row_sector':'sector'},inplace=True)
trd = pd.merge(left=trd,right=cons,how='left',on=['region','sector'])

# merge on gdp
gdp = output.groupby(['col_region'])['VA'].sum().reset_index()
gdp.rename(columns={'col_region':'region','VA':'GDP'},inplace=True)
trd = pd.merge(left=trd,right=gdp,how='left',on=['region'])
trd.loc[trd.sector=='TOT','VA']=trd.loc[trd.sector=='TOT','GDP']

# restrict data to USA
df=trd.loc[(trd.sector!='TOT') & (trd.region=='1-USA'),:].reset_index(drop=True)
df['trd']=df.ex+df.im
df['trd_M']=df.ex_M+df.im_M
df['nx'] = df.ex-df.im
df['nx_M'] = df.ex_M-df.im_M
df['nx_F'] = df.ex_F-df.im_F
df['nx_C'] = df.ex_C-df.im_C
df['nx_I'] = df.ex_I-df.im_I
df['trd2'] = df.trd
df['im2'] = df.im
df['ex2'] = df.ex
df['nx2'] = df.nx
df['im2_M'] = df.im_M
df['ex2_M'] = df.ex_M
df['nx2_M'] = df.nx_M
df['im2_F'] = df.im_F
df['ex2_F'] = df.ex_F
df['nx2_F'] = df.nx_F
df['im2_C'] = df.im_C
df['ex2_C'] = df.ex_C
df['nx2_C'] = df.nx_C
df['im2_I'] = df.im_I
df['ex2_I'] = df.ex_I
df['nx2_I'] = df.nx_I


# normalize by sectoral value added and GDP
for col in df.columns:
    if col in ['ex_M','ex_F','ex_C','ex_I','ex',
               'im_M','im_F','im_C','im_I','im',
               'nx_M','nx_F','nx_C','nx_I','nx']:
        df[col]=100*df[col]/df.GO
    elif col in ['ex2_M','ex2_F','ex2_C','ex2_I','ex2',
                 'im2_M','im2_F','im2_C','im2_I','im2',
                 'nx2_M','nx2_F','nx2_C','nx2_I','nx2']:
        df[col]=100*df[col]/df.GDP

sectors = {'1-UPSTREAM-HI':'Oil',
           '2-UPSTREAM-LO':'Steel',
           '3-DOWNSTREAM-HI':'Toys',
           '4-DOWNSTREAM-LO':'Cars',
           '5-SERVICES':'Svcs',
           '6-CONSTRUCTION':'Const'}

slist=list(sectors.keys())
snames=[sectors[s] for s in sectors.keys()]

partners = {'2-CHN':'China',
            '3-ROW':'Rest of world'}

pnames = [partners[p] for p in partners.keys()]

#-------------------------------------------------------------------
# stacked bar plots of trade composition across countries, by sector

cols=['ex','im','nx','ex2','im2','nx2']
ylims=[(0,80),(0,80),(-80,10),(0,6),(0,6),(-2,1.5)]
df2 = df[df.partner!='TOT'].groupby(['sector','partner'])[cols].sum().reset_index()
width=0.75

fcnt=0
for c in cols:
    fig,ax=create_fig()

    data = {}
    for partner,name in partners.items():
        data[name] = np.array([df2[c][(df2.partner==partner) & (df2.sector==s)].values[0] for s in slist[0:-1]])

    bottom_pos=np.zeros(len(slist)-1)
    bottom_neg=np.zeros(len(slist)-1)

    bcnt=0
    for name, values in data.items():
        bottom=which_bottom(bottom_pos,bottom_neg,values)
        ax.bar(snames[:-1], values, width, label=name, bottom=bottom, color=colors[bcnt], alpha=0.99, edgecolor='black', lw=1)
        bottom_pos=bottom_pos+vmax(values)
        bottom_neg=bottom_neg+vmin(values)
        bcnt += 1

    if(bottom_neg.min()<1e-6):
        ax.axhline(0.0,color='black',ls='-',lw=1)
               
    ax.set_ylim(ylims[fcnt][0],ylims[fcnt][1])
    
    #    if fcnt in [1,4]:
    ax.legend(frameon=True, loc='best',fontsize=8, framealpha=1, edgecolor='white',borderpad=0,borderaxespad=1)
        
    fig.tight_layout()    
    plt.savefig('output/fig_data_sectoral_trade_by_region_%s.pdf'%c)
    
    plt.clf()
    fcnt +=1

plt.close('all')

#-------------------------------------------------------------------
# stacked bar plots of trade composition across uses, by sector

cols=['ex_M','im_M','nx_M','ex_F','im_F','nx_F','ex_C','im_C','nx_C','ex_I','im_I','nx_I',
      'ex2_M','im2_M','nx2_M','ex2_F','im2_F','nx2_F','ex2_C','im2_C','nx2_C','ex2_I','im2_I','nx2_I']

ylims=[(0,80),(0,80),(-80,10),(0,6),(0,6),(-2,1.5)]
df2 = df[df.partner!='TOT'].groupby(['sector'])[cols].sum().reset_index()
width=0.75

#uses={'Intermediate':'_M','Final':'_F'}
uses={'Intermediate':'_M','Investment':'_I','Consumption':'_C'}

fcnt=0
for c in ['ex','im','nx','ex2','im2','nx2']:
    fig,ax=create_fig()

    data = {}
    for name,suff in uses.items():
        data[name] = np.array([df2[c+suff][(df2.sector==s)].values[0] for s in slist[0:-1]])

    bottom_pos=np.zeros(len(slist)-1)
    bottom_neg=np.zeros(len(slist)-1)

    bcnt=0
    for name, values in data.items():
        bottom=which_bottom(bottom_pos,bottom_neg,values)
        ax.bar(snames[:-1], values, width, label=name, bottom=bottom, color=colors[bcnt], alpha=0.99, edgecolor='black', lw=1)
        bottom_pos=bottom_pos+vmax(values)
        bottom_neg=bottom_neg+vmin(values)
        bcnt += 1

    if(bottom_neg.min()<1e-6):
        ax.axhline(0.0,color='black',ls='-',lw=1)
               
    ax.set_ylim(ylims[fcnt][0],ylims[fcnt][1])

    #if fcnt in [1,4]:
    ax.legend(frameon=True, loc='best',fontsize=8, framealpha=1, edgecolor='white',borderpad=0,borderaxespad=1)
        
    fig.tight_layout()    
    plt.savefig('output/fig_data_sectoral_trade_by_use_%s.pdf'%c)
    
    plt.clf()
    fcnt +=1

plt.close('all')


#---------------------------------------------------------------------
# Downstream linkages (direct requirement coefficients)
    
M = intermediates.groupby(['col_region','col_sector','row_region','row_sector'])['M'].sum().reset_index()
M = pd.merge(left=M,right=output,how='left',on=['col_region','col_sector'])
M['DR'] = 100*M.M/M.GO

Mt = M.loc[(M.col_region=='1-USA')].groupby(['col_sector','row_sector'])['DR'].sum().reset_index()
Md = M.loc[(M.col_region=='1-USA') & (M.col_region==M.row_region)].groupby(['col_sector','row_sector'])['DR'].sum().reset_index()
Mf = M.loc[(M.col_region=='1-USA') & (M.col_region!=M.row_region)].groupby(['col_sector','row_sector'])['DR'].sum().reset_index()
Mc = M.loc[(M.col_region=='1-USA') & (M.row_region=='2-CHN')].groupby(['col_sector','row_sector'])['DR'].sum().reset_index()
Mr = M.loc[(M.col_region=='1-USA') & (M.row_region=='3-ROW')].groupby(['col_sector','row_sector'])['DR'].sum().reset_index()

M = [Mt,Md,Mf,Mc,Mr]
for i in range(5):

    df = M[i]
    fig,ax=create_fig()
    
    data = {}
    for s_src, sn in sectors.items():
        if(sn!='Svcs' and sn!='Const'):
            data[sn] = np.array( [df.DR[(df.col_sector==s_dest) & (df.row_sector==s_src)].values[0] for s_dest in slist] )
        
    bottom_pos=np.zeros(len(slist))
    bottom_neg=np.zeros(len(slist))
    
    bcnt=0
    for name, values in data.items():
        bottom=which_bottom(bottom_pos,bottom_neg,values)
        ax.bar(snames, values, width, label=name, bottom=bottom, color=colors[bcnt], alpha=0.99, edgecolor='black', lw=1)
        bottom_pos=bottom_pos+vmax(values)
        bottom_neg=bottom_neg+vmin(values)
        bcnt += 1

    if(bottom_neg.min()<1e-6):
        ax.axhline(0.0,color='black',ls='-',lw=1)
               
    ax.set_ylim(0,40)

    if(i==0):
        ax.legend(frameon=True, loc='upper right',fontsize=8, framealpha=1, edgecolor='white',borderpad=0,borderaxespad=1)
        
    fig.tight_layout()
    plt.savefig('output/fig_data_linkages_downstream_%d.pdf'%i)
    
    plt.clf()
    fcnt +=1


# -------------------------------------------------------------------------------------------------------
# Upstream linkages

M = intermediates.groupby(['col_region','col_sector','row_region','row_sector'])['M'].sum().reset_index()
M = pd.merge(left=M,right=output.rename(columns={'col_region':'row_region','col_sector':'row_sector'}),
             how='left',on=['row_region','row_sector'])
M['UR'] = 100*M.M/M.GO

Mt = M.loc[(M.col_region=='1-USA')].groupby(['col_sector','row_sector'])['UR'].sum().reset_index()
Md = M.loc[(M.col_region=='1-USA') & (M.col_region==M.row_region)].groupby(['col_sector','row_sector'])['UR'].sum().reset_index()
Mf = M.loc[(M.col_region=='1-USA') & (M.col_region!=M.row_region)].groupby(['col_sector','row_sector'])['UR'].sum().reset_index()
Mc = M.loc[(M.col_region=='1-USA') & (M.row_region=='2-CHN')].groupby(['col_sector','row_sector'])['UR'].sum().reset_index()
Mr = M.loc[(M.col_region=='1-USA') & (M.row_region=='3-ROW')].groupby(['col_sector','row_sector'])['UR'].sum().reset_index()

M = [Mt,Md,Mf,Mc,Mr]
for i in range(5):

    df = M[i]
    fig,ax=create_fig()
    
    data = {}
    for s_dest, sn in sectors.items():
        #if(sn!='Services' and sn!='Construction'):
        data[sn] = np.array( [df.UR[(df.col_sector==s_dest) & (df.row_sector==s_src)].values[0] for s_src in slist[0:-2]] )
        
    bottom_pos=np.zeros(len(slist[0:-2]))
    bottom_neg=np.zeros(len(slist[0:-2]))
    
    bcnt=0
    for name, values in data.items():
        bottom=which_bottom(bottom_pos,bottom_neg,values)
        ax.bar(snames[0:-2], values, width, label=name, bottom=bottom, color=colors[bcnt], alpha=0.99, edgecolor='black', lw=1)
        bottom_pos=bottom_pos+vmax(values)
        bottom_neg=bottom_neg+vmin(values)
        bcnt += 1

    if(bottom_neg.min()<1e-6):
        ax.axhline(0.0,color='black',ls='-',lw=1)
               
    ax.set_ylim(0,80)

    if(i==0):
        ax.legend(frameon=True, loc='upper right',fontsize=8, framealpha=1, edgecolor='white',borderpad=0,borderaxespad=1)
        
    fig.tight_layout()
    plt.savefig('output/fig_data_linkages_upstream_%d.pdf'%i)
    
    plt.clf()
    fcnt +=1


