import numpy as np
import pandas as pd
import matplotlib as mpl
import pandas as pd
import locale
locale.setlocale(locale.LC_ALL,'en_US.utf8')

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


# aggregate within services sector
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
U = U.loc[~(U.row_sector.isin(['SVCS','CONS'])),:]

# display results
U = U.rename(columns={'row_sector':'ind'})
names = pd.read_csv('../../data/industry_names.csv')
U = pd.merge(left=U,right=names,how='left',on='ind')
U2 = U.groupby(['ind','name'])['U'].mean().sort_values(ascending=False).reset_index()
U2['upstream'] = np.where(U2.U > U2.U.median(),True,False)
print(U2)

#############################################################################
print('\nAggregating across regions and sectors...')

# link back to main dataframe and finish aggregating across sectors
agged = pd.merge(left=agged,
                 right=U2.rename(columns={'ind':'row_sector','upstream':'row_upstream'})[['row_sector','row_upstream']],
                 how='left',on='row_sector')

agged = pd.merge(left=agged,
                 right=U2.rename(columns={'ind':'col_sector','upstream':'col_upstream'})[['col_sector','col_upstream']],
                 how='left',on='col_sector')
agged.loc[agged.col_upstream==True,'col_sector']='UP'
agged.loc[agged.col_upstream==False,'col_sector']='DN'
agged.loc[agged.row_upstream==True,'row_sector']='UP'
agged.loc[agged.row_upstream==False,'row_sector']='DN'
agged = agged.drop(['col_upstream','row_upstream'],axis=1)
agged = agged.groupby(['row_region','row_sector','col_use','col_region','col_sector'])['value'].sum().reset_index()
agged = agged.sort_values(by=['row_region','row_sector','col_use','col_region','col_sector'])

# order region and sectors as desired
agged.loc[agged.row_region=='USA','row_region'] = '1-USA'
agged.loc[agged.row_region=='CHN','row_region'] = '2-CHN'
agged.loc[agged.row_region=='ROW','row_region'] = '3-ROW'
agged.loc[agged.col_region=='USA','col_region'] = '1-USA'
agged.loc[agged.col_region=='CHN','col_region'] = '2-CHN'
agged.loc[agged.col_region=='ROW','col_region'] = '3-ROW'

agged.loc[agged.row_sector=='UP','row_sector'] = '1-UPSTREAM'
agged.loc[agged.row_sector=='DN','row_sector'] = '2-DOWNSTREAM'
agged.loc[agged.row_sector=='SVCS','row_sector'] = '3-SERVICES'
agged.loc[agged.row_sector=='CONS','row_sector'] = '4-CONSTRUCTION'
agged.loc[agged.col_sector=='UP','col_sector'] = '1-UPSTREAM'
agged.loc[agged.col_sector=='DN','col_sector'] = '2-DOWNSTREAM'
agged.loc[agged.col_sector=='SVCS','col_sector'] = '3-SERVICES'
agged.loc[agged.col_sector=='CONS','col_sector'] = '4-CONSTRUCTION'

# separate into main components in same structure as in NAFTA paper
intermediates = agged[(agged.col_region!='TOT') &
                      (agged.row_region!='TOT') &
                      (agged.col_use=='INT')].groupby(['col_region','col_sector','row_region','row_sector'])['value']\
                      .sum().reset_index().rename(columns={'value':'M'})

consumption = agged[(agged.col_sector=='CC') &
                    (agged.col_region!='TOT') &
                    (agged.row_region!='TOT')].groupby(['col_region','row_region','row_sector'])['value']\
                    .sum().reset_index().rename(columns={'value':'C'})

investment = agged[(agged.col_sector=='INV') &
                    (agged.col_region!='TOT') &
                    (agged.row_region!='TOT')].groupby(['col_region','row_region','row_sector'])['value']\
                    .sum().reset_index().rename(columns={'value':'I'})

value_added = agged[(agged.col_region!='TOT') &
                    (agged.col_sector!='TOT') &
                    (agged.col_use=='INT') &
                    (agged.row_sector=='VA')].groupby(['col_region','col_sector'])['value']\
                    .sum().reset_index().rename(columns={'value':'VA'})

gross_output = agged[(agged.col_region!='TOT') &
                     (agged.col_sector!='TOT') &
                     (agged.col_use=='INT') &
                     (agged.row_sector=='GO')].groupby(['col_region','col_sector'])['value']\
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
final_demand.loc[ (final_demand.row_sector=='4-CONSTRUCTION') & (final_demand.col_region != final_demand.row_region), 'I']=0
intermediates.loc[ (intermediates.row_sector=='4-CONSTRUCTION') , 'M']=0
final_demand.loc[ (final_demand.row_sector=='4-CONSTRUCTION'), 'C']=0


    
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
sectors = ['$T_U$','$T_D$','$T_S$','$N_C$']

def write_iomat_csv(iomat,fname):
    usgdp = iomat[-1,0:ns].sum()
    iomat2 = np.vstack((iomat,np.sum(iomat,axis=0).reshape((1,nc*ns+nc*2))))
    iomat2 = np.hstack((iomat2,np.sum(iomat2,axis=1).reshape((nc*ns+2,1))))
    iomat2 = 100*iomat2/usgdp
    np.savetxt(fname=fname+'.csv',X=iomat2,fmt='%0.15f',delimiter=' ')

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
        #file.write('\\begin{landscape}\n')
        #file.write('\\begin{table}[p]\n')
        #file.write('\\renewcommand{\\arraystretch}{1.2}\n')
        #file.write('\\begin{center}\n')
        #file.write('\\caption{'+caption+', intermediate inputs portion '+units+'}\n')
        #file.write('\\label{tab:'+label+'_m}\n')
        #file.write('\\footnotesize\n')
        
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
                    if Fx[i*ns+ii,j]>1e-6:
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
va = output.groupby(['col_region','col_sector'])['VA'].sum().reset_index()
va.rename(columns={'col_region':'region','col_sector':'sector'},inplace=True)
trd = pd.merge(left=trd,right=va,how='left',on=['region','sector'])

# merge on consumption
cons = final_demand.groupby(['col_region','row_sector'])['C'].sum().reset_index()
cons.rename(columns={'col_region':'region','row_sector':'sector'},inplace=True)
trd = pd.merge(left=trd,right=cons,how='left',on=['region','sector'])

# merge on gdp
gdp = output.groupby(['col_region'])['VA'].sum().reset_index()
gdp.rename(columns={'col_region':'region','VA':'GDP'},inplace=True)
trd = pd.merge(left=trd,right=gdp,how='left',on=['region'])
trd.loc[trd.sector=='TOT','VA']=trd.loc[trd.sector=='TOT','GDP']

# make latex table
trd=trd.groupby(['region','partner','sector']).mean().reset_index()

sector_names = {'1-UPSTREAM':'Upstream goods',
                '2-DOWNSTREAM':'Downstream goods',
                '3-SERVICES':'Services',
                '4-CONSTRUCTION':'Construction',
                'TOT':'Total'}

country_names = {'1-USA':'United States',
                 '2-CHN':'China',
                 '3-ROW':'Rest of world',
                 'TOT':'Total'}
partners = {'1-USA':['TOT','2-CHN','3-ROW'],'2-CHN':['TOT','1-USA','3-ROW'],'3-ROW':['TOT','1-USA','2-CHN']}
panels = {'1-USA':'(a)','2-CHN':'(b)','3-ROW':'(c)'}

with open('output/icio_summary.tex','w') as file:
    #file.write('\\begin{table}[p]\n')
    #file.write('\\renewcommand{\\arraystretch}{1.2}\n')
    #file.write('\\begin{center}\n')
    #file.write("\\caption{Sectoral production and trade in NAFTA (2014 data, percent GDP)}\n")
    #file.write('\\label{tab:key_facts}\n')
    #file.write('\\footnotesize\n')
    file.write('\\begin{tabular}{lccccc}\n')
    file.write('\\toprule\n')
    file.write('\\multicolumn{1}{p{2cm}}{\\centering Quantity} & ')

    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Upstream\\\\goods} & ')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Downstream\\\\goods} &')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Services} &')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Construction} & ')
    file.write('\\multicolumn{1}{p{1.5cm}}{\\centering Total}\\\\\n')
    file.write('\\midrule\n')

    for c in ['1-USA','2-CHN','3-ROW']:
        file.write('\\multicolumn{5}{l}{\\textit{'+panels[c]+' '+country_names[c]+'}}\\\\\n')
        mask=trd.region==c

        file.write('Value added')
        mask2 = np.logical_and(mask,trd.partner=='TOT')
        for s in sector_names.keys():
            mask3=np.logical_and(mask2,trd.sector==s)
            masked=trd[mask3]
            val = 100.0*masked['VA']/masked['GDP']
            file.write('& %0.2f' % val.iloc[0])
        file.write('\\\\\n')
        
        for p in partners[c]:
            mask2=np.logical_and(mask,trd.partner==p)
            if p=='TOT':
                file.write('Exports')
            else:
                file.write('\quad to ' + country_names[p])
            for s in sector_names.keys():
                mask3=np.logical_and(mask2,trd.sector==s)
                masked=trd[mask3]
                val = 100.0*(masked['ex'])/masked['GDP']
                file.write('& %0.2f' % val.iloc[0])            
            file.write('\\\\\n')

        for p in partners[c]:
            mask2=np.logical_and(mask,trd.partner==p)
            if p=='TOT':
                file.write('Imports')
            else:
                file.write('\quad from ' + country_names[p])
            for s in sector_names.keys():
                mask3=np.logical_and(mask2,trd.sector==s)
                masked=trd[mask3]
                val = 100.0*(masked['im'])/masked['GDP']
                file.write('& %0.2f' % val.iloc[0])            
            file.write('\\\\\n')

        for p in partners[c]:
            mask2=np.logical_and(mask,trd.partner==p)
            if p=='TOT':
                file.write('Net exports')
            else:
                file.write('\quad with ' + country_names[p])
            for s in sector_names.keys():
                mask3=np.logical_and(mask2,trd.sector==s)
                masked=trd[mask3]
                val = 100.0*masked['tb']/masked['GDP']
                file.write('& %0.2f' % val.iloc[0])            
            file.write('\\\\\n')

        if(c!='ROW'):
            file.write('\\\\\n')

    file.write('\\bottomrule\n')
    file.write('\\end{tabular}\n')
#    file.write('\\normalsize\n')
#    file.write('\\end{center}\n')
#    file.write('\\end{table}\n')
