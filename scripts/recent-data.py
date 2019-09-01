#####################################################################################
# imports, small functions, etc.

import numpy as np
from scipy import signal
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':26})
mpl.rc('font',size=26)
mpl.rc('text', usetex=True)
mpl.rc('lines',linewidth=1.5)
#mpl.rc('savefig',transparent=False)
mpl.rc('savefig',bbox='tight')
mpl.rc('savefig',format='pdf')

import locale
locale.setlocale(locale.LC_ALL,'en_US.utf8')

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':26})
mpl.rc('font',size=26)
mpl.rc('text', usetex=True)
mpl.rc('lines',linewidth=1.5)
#mpl.rc('savefig',transparent=False)
mpl.rc('savefig',bbox='tight')
mpl.rc('savefig',format='pdf')

def yrmo_to_yrq(yrmo):
    yr=yrmo[0:4]
    m=int(yrmo[5:7])
    q=0
    if m<=3:
        q=1
    elif m<=6:
        q=2
    elif m<=9:
        q=3
    else:
        q=4
    return yr+'Q'+str(q)

def quarterize(df,cols,sum=True):
    df['Period']=df.Month.apply(yrmo_to_yrq)

    if sum==True:
        return df.groupby('Period')[cols].sum().reset_index()
    else:
        return df.groupby('Period')[cols].mean().reset_index()

def log(df,cols):
    for col in cols:
        df[col+'_l'] = np.log(df[col]/df[col][0])

def detrend(df,cols,endpoint):
    for col in cols:
        g = (df[col][endpoint]/df[col][0])**(1.0/len(df[col][0:endpoint]))
        df['g'] = 1.0
        df.loc[1:,'g'] = g
        df['cg'] = df['g'].cumprod()
        df[col+'_t'] = df[col][0]*df['cg']
        df[col+'_dt'] = 100.0*(df[col]/df[col+'_t']-1.0)

def ticks(ax,periods,step):
    n = len(periods)
    tt=range(0,n,step)
    mask=[i in tt for i in periods.index.values]
    labels=periods[mask].reset_index(drop=True)
    ax.set_xlim(periods.index.values[0],periods.index.values[-1])
    ax.set_xticks(periods.index.values,minor=True)
    ax.set_xticks(tt,minor=False)
    ax.set_xticklabels(labels,minor=False,rotation=45)
    

#####################################################################################

# national accounts
rna = pd.read_csv('../data/eurostat-national-accounts-real.csv')
nna = pd.read_csv('../data/eurostat-national-accounts-nominal.csv')

# bilateral trade data
eutrd = pd.read_csv('../data/eurostat-eu-trd-ner.csv')
eutrd['EU_exports_gbp']=eutrd.EU_exports_euro * eutrd.gbp_euro
eutrd['EU_imports_gbp']=eutrd.EU_imports_euro * eutrd.gbp_euro
eutrd_q = quarterize(eutrd,['EU_exports_gbp','EU_imports_gbp','EU_exports_euro','EU_imports_euro'],sum=True)
eutrd_q=eutrd_q[eutrd_q.Period!='2018Q2']

# eu nominal gdp
eungdp = pd.read_csv('../data/eu27_ngdp.csv')

# RER
wrer = pd.read_csv('../data/Weighted_RER.csv')
wrer_eu_q = quarterize(wrer.loc[wrer.Partner==1],['RER'],sum=False)
wrer_rw_q = quarterize(wrer.loc[wrer.Partner==2],['RER'],sum=False)
wrer_eu_q['RER'] = wrer_eu_q.RER/wrer_eu_q.RER[wrer_eu_q.Period=='2015Q2'].values[0]
wrer_rw_q['RER'] = wrer_rw_q.RER/wrer_rw_q.RER[wrer_rw_q.Period=='2015Q2'].values[0]

# merge with other nominal data
nna = pd.merge(left=nna,right=eutrd_q,how='left',on='Period')
nna = pd.merge(left=nna,right=eungdp,how='left',on='Period')

# merge real data with nominal variables and normalize by GDP
nna['TB_pctY'] = 100.0*(nna['EX']-nna['IM'])/nna['GDP']
nna['EX_pctY'] = 100.0*(nna['EX'])/nna['GDP']
nna['IM_pctY'] = 100.0*(nna['IM'])/nna['GDP']
nna['EX_EU_pctY'] = 100.0*(nna['EU_exports_gbp'])/nna['GDP']
nna['IM_EU_pctY'] = 100.0*(nna['EU_imports_gbp'])/nna['GDP']
nna['EX_EU_pctY_euro'] = 100.0*(nna['EU_exports_euro'])/nna['EU_GDP_euro']
nna['IM_EU_pctY_euro'] = 100.0*(nna['EU_imports_euro'])/nna['EU_GDP_euro']
nna['GFCF_pctY'] = 100.0*(nna['GFCF'])/nna['GDP']
rna = pd.merge(left=rna,
               right=nna[['Period','GFCF_pctY','EX_pctY','IM_pctY','EX_EU_pctY','IM_EU_pctY','EX_EU_pctY_euro','IM_EU_pctY_euro','TB_pctY']],
               how='left',on='Period')

rna = rna[148:].reset_index()

file=open('tex/recent-data.tex','wb')

file.write('\\begin{table}[h!]\n')
file.write('\\renewcommand{\\arraystretch}{1.2}\n')
file.write('\\begin{center}\n')

file.write('\\caption{Recent U.K. macroeconomic and trade dynamics}\n')
file.write('\\label{tab:recent-data}\n')

# number of columns
file.write('\\begin{tabular}{l c c c}\n')
file.write('\\toprule\n')

# headers
file.write('Variable & 2012Q1--2015Q2 & 2015Q3--2016Q2 & 2016Q3--2018Q2\\\\\n')
file.write('\\midrule\n')

# data
def growth(col):

    tmpdbl = 100.*((rna[col][13]/rna[col][0])**(4./13.) - 1.0)
    tmpstr = locale.format('%0.2f',tmpdbl,grouping=True)
    file.write('&'+tmpstr)

    tmpdbl = 100.*((rna[col][17]/rna[col][13])**(4./4.) - 1.0)
    tmpstr = locale.format('%0.2f',tmpdbl,grouping=True)
    file.write('&'+tmpstr)

    tmpdbl = 100.*((rna[col][25]/rna[col][17])**(4./8.) - 1.0)
    tmpstr = locale.format('%0.2f',tmpdbl,grouping=True)
    file.write('&'+tmpstr+'\\\\\n')

def avg(col):

    tmpdbl = rna[col][0:14].mean()
    tmpstr = locale.format('%0.2f',tmpdbl,grouping=True)
    file.write('&'+tmpstr)

    tmpdbl = rna[col][14:18].mean()
    tmpstr = locale.format('%0.2f',tmpdbl,grouping=True)
    file.write('&'+tmpstr)

    tmpdbl = rna[col][18:26].mean()
    tmpstr = locale.format('%0.2f',tmpdbl,grouping=True)
    file.write('&'+tmpstr+'\\\\\n')

file.write('\\multicolumn{4}{l}{\\textit{(a) National income accounts}}\\\\\n')
file.write('Real GDP growth (pct. per year)')
growth('GDP')

#file.write('Employment growth (pct. per year)')
#growth('EMP')

file.write('Consumption growth (pct. per year)')
growth('C')

file.write('Investment (pct. GDP)')
avg('GFCF_pctY')

file.write('Net exports (pct. GDP)')
avg('TB_pctY')

file.write('\\\\\n\\multicolumn{4}{l}{\\textit{(b) International trade}}\\\\\n')
file.write('Exports (pct. GDP)')
avg('EX_pctY')

file.write('Imports (pct. GDP)')
avg('IM_pctY') 

file.write('Goods exports to E.U. (pct. GDP)')
avg('EX_EU_pctY')

file.write('Goods imports from E.U. (pct. GDP)')
avg('IM_EU_pctY') 

file.write('\\bottomrule\n')
file.write('\\end{tabular}\n')
file.write('\\end{center}\n')
file.write('\\end{table}\n')
file.close()

#####################################################################################
# plot

rna=rna[0:-1].reset_index()

colors=['#377eb8','#4daf4a','#e41a1c','#1b9e77','#d95f02']

fig=plt.figure()
ax=fig.add_axes([0,0,1,1])       

ax.plot(np.log(rna.GDP/rna.GDP[0]),linewidth=3,alpha=0.75,color=colors[0])
ax.plot(np.log(rna.C/rna.C[0]),linewidth=3,alpha=0.75,color=colors[1],dashes=(3,3))
ax.plot(np.log(rna.GFCF/rna.GFCF[0]),linewidth=3,alpha=0.75,color=colors[2],dashes=(9,6))
ax.axvline(13,color='black',linewidth=1)
ax.axvline(17,color='black',linewidth=1)
ax.set_xticks(range(len(rna.Period)))
ticklabs=['' for x in ax.get_xticks()]
ticklabs[0]='2012Q1'
ticklabs[13]='2015Q2'
ticklabs[17]='2016Q2'
ax.set_xticklabels(ticklabs)
ax.annotate('Referendum announced',
            size=20,
            xy=(12.9,0.16),
            xytext=(2,.17),
            arrowprops=dict(width=0.5,headwidth=5,frac=0.2))
ax.annotate('Brexit vote',
            size=20,
            xy=(17.1,-0.02),
            xytext=(18,-0.02),
            arrowprops=dict(width=0.5,headwidth=5,frac=0.2))
ax.annotate('GDP',size=20,xy=(22,0.1225),xytext=(22,0.1225))
ax.annotate('Consumption',size=20,xy=(18.5,0.0775),xytext=(18.5,0.0775))
ax.annotate('Investment',size=20,xy=(4.5,-0.03),xytext=(4.5,-0.03))
#plt.ylim(-0.1,0.2)
#plt.legend(['GDP','Consumption','Investment'],loc='lower center',prop={'size':20},ncol=3)
plt.xlim(0,len(rna.Period)-1)
plt.savefig('fig/recent-data-na.pdf')
ax.set_ylabel('log change (base period: 2012Q1)')
ax.set_title("UK nat'l accts since Brexit referendum",y=1.05)
plt.savefig('fig/recent-data-na-web.png')
plt.clf()
plt.close('all')

fig=plt.figure()
ax=fig.add_axes([0,0,1,1])       

ax.plot(rna.TB_pctY,linewidth=3,alpha=0.75,color=colors[0])
ax.plot(rna.EX_pctY,linewidth=3,alpha=0.75,color=colors[1],dashes=(3,3))
ax.plot(rna.IM_pctY,linewidth=3,alpha=0.75,color=colors[2],dashes=(9,6))
ax.plot(rna.EX_EU_pctY,linewidth=3,alpha=0.75,color=colors[3],dashes=(6,3,6,3))
ax.plot(rna.IM_EU_pctY,linewidth=3,alpha=0.75,color=colors[4],dashes=(3,1,3,1))
ax.axvline(13,color='black',linewidth=1)
ax.axvline(17,color='black',linewidth=1)
ax.set_xticks(range(len(rna.Period)))
ticklabs=['' for x in ax.get_xticks()]
ticklabs[0]='2012Q1'
ticklabs[13]='2015Q2'
ticklabs[17]='2016Q2'
ax.set_xticklabels(ticklabs)
#ax.annotate('Referendum announced',
#            size=20,
#            xy=(12.9,20),
#            xytext=(2,20),
#            arrowprops=dict(width=0.5,headwidth=5,frac=0.2))
#ax.annotate('Brexit vote',
#            size=20,
#            xy=(17.1,20),
#            xytext=(18,20),
#            arrowprops=dict(width=0.5,headwidth=5,frac=0.2))
ax.annotate('Imports',size=20,xy=(4,32.5),xytext=(4,32.5))
ax.annotate('Exports',size=20,xy=(1.5,27.5),xytext=(1.5,27.5))
ax.annotate('Net exports',size=20,xy=(19,0),xytext=(19,0))
ax.annotate('EU goods imports',size=20,xy=(3,12.5),xytext=(3,12.5))
ax.annotate('EU goods exports',size=20,xy=(1,5.0),xytext=(1,5.0))
#plt.legend(['Net exports','Exports','Imports','EU goods exports','EU goods imports'],loc='lower center',prop={'size':20},ncol=2)
plt.ylim(-5,35)
plt.xlim(0,len(rna.Period)-1)
plt.savefig("fig/recent-data-trd.pdf")
ax.set_ylabel('pct. UK GDP')
ax.set_title("UK trade flows since Brexit referendum",y=1.05)
plt.savefig('fig/recent-data-trd-web.png')
plt.clf()

fig=plt.figure()
ax=fig.add_axes([0,0,1,1])       
ax.plot(rna.EX_EU_pctY_euro,linewidth=3,alpha=0.75,color=colors[0])
ax.plot(rna.IM_EU_pctY_euro,linewidth=3,alpha=0.75,color=colors[1],dashes=(9,6))
ax.axvline(13,color='black',linewidth=1)
ax.axvline(17,color='black',linewidth=1)
ax.set_xticks(range(len(rna.Period)))
ticklabs=['' for x in ax.get_xticks()]
ticklabs[0]='2012Q1'
ticklabs[13]='2015Q2'
ticklabs[17]='2016Q2'
ax.set_xticklabels(ticklabs)
#ax.annotate('Referendum announced',
#            size=20,
#            xy=(12.9,2),
#            xytext=(2,2),
#            arrowprops=dict(width=0.5,headwidth=5,frac=0.2))
#ax.annotate('Brexit vote',
#            size=20,
#            xy=(17.1,2),
#            xytext=(18,2),
#            arrowprops=dict(width=0.5,headwidth=5,frac=0.2))
ax.annotate('EU goods imports',size=20,xy=(3,2.5),xytext=(3,2.5))
ax.annotate('EU goods exports',size=20,xy=(1,1.5),xytext=(1,1.5))
plt.xlim(0,len(rna.Period)-1)
plt.savefig('fig/recent-data-trd-euro.pdf')
plt.clf()

fig=plt.figure()
ax=fig.add_axes([0,0,1,1])       
ax.plot(wrer_eu_q.RER[4:-1],linewidth=3,alpha=0.75,color=colors[0])
ax.plot(wrer_rw_q.RER[4:-1],linewidth=3,alpha=0.75,color=colors[1],dashes=(9,6))
ax.axvline(13,color='black',linewidth=1)
ax.axvline(17,color='black',linewidth=1)
ax.set_xticks(range(len(rna.Period)))
ticklabs=['' for x in ax.get_xticks()]
ticklabs[0]='2012Q1'
ticklabs[13]='2015Q2'
ticklabs[17]='2016Q2'
ax.set_xticklabels(ticklabs)
#ax.annotate('Referendum announced',
#            size=20,
#            xy=(12.9,0.5),
#            xytext=(2,0.5),
#            arrowprops=dict(width=0.5,headwidth=5,frac=0.2))
#ax.annotate('Brexit vote',
#            size=20,
#            xy=(17.1,0.5),
#            xytext=(18,0.5),
#            arrowprops=dict(width=0.5,headwidth=5,frac=0.2))
ax.annotate('E.U. RER',size=20,xy=(4,1.2),xytext=(4,1.2))
ax.annotate('R.W. RER',size=20,xy=(8,0.925),xytext=(8,0.925))
plt.xlim(0,len(rna.Period)-1)
plt.savefig('fig/recent-data-rer.pdf')
plt.clf()




plt.close('all')

