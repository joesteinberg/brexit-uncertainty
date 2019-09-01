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
mpl.rc('savefig',bbox='tight')
mpl.rc('savefig',format='pdf')  

#####################################################################################

colors=['#377eb8','#4daf4a','#e41a1c','#1b9e77','#d95f02']

trd = pd.read_csv('../data/comtrade-uk-olddata-efta.csv')
gdp = pd.read_csv('../data/fred_ukgdp_usd.csv')
merged = pd.merge(left=trd,right=gdp,how='left',on='Year')

wld=trd[trd['Partner ISO']=='WLD'].drop('Partner ISO',axis=1)
wld.rename(columns={'Trade Value (US$)':'TOT'},inplace=True)
merged = pd.merge(left=merged,right=wld,how='left',on=['Year','Trade Flow'])

merged['frac'] = 100*merged['Trade Value (US$)']/merged.UK_GDP_USD
merged['frac2'] = 100*merged['Trade Value (US$)']/merged.TOT

efta = merged[merged['Partner ISO']!='WLD'].groupby(['Year','Trade Flow'])[['frac','frac2']].sum().reset_index()

fig=plt.figure()
ax=fig.add_axes([0,0,1,1])       

yrs=efta.Year.unique()

y1=efta.loc[efta['Trade Flow']=='Export','frac2'].values
y2=efta.loc[efta['Trade Flow']=='Import','frac2'].values
ax.plot(yrs,y1+y2,linewidth=3,alpha=0.75,color=colors[0])
ax.axvline(1973,color='black',linewidth=1)
ax.set_xlim(1962,1992)
ax.annotate('U.K. leaves EFTA',
            size=20,
            xy=(1973.25,24),
            xytext=(1975,25),
            arrowprops=dict(width=0.5,headwidth=5,frac=0.2))
plt.savefig('fig/efta_trd.pdf')

plt.clf()


##########################################################################3

print 'Drop in EFTA share of UK trade: % 0.2f p.p., or %0.2f pct' % ( ((min(y1+y2))-(max(y1+y2))),
                                                                      100*((min(y1+y2))/(max(y1+y2))-1.0) )

path = '../quanal/dyn_mkt_pen/output/'
models=[]
models.append(pd.read_csv(path+'vars_nobrexit_uk_baseline.csv'))
models.append(pd.read_csv(path+'vars_stoch_soft_uk_baseline.csv'))
models.append(pd.read_csv(path+'vars_stoch_hard_uk_baseline.csv'))

for m in models:
    m['tot_trd']=m.ex1+m.im1+m.ex2+m.im2
    m['eu_share']=(m.ex1+m.im1)/m.tot_trd

print 'Drop in EU share of UK trade (soft): %0.2f p.p., or %0.2f pct' % ( 100*(models[1].eu_share[49]-models[1].eu_share[0]),
                                                                          100*(models[1].eu_share[49]/models[1].eu_share[0]-1.0) )
print 'Drop in EU share of UK trade (hard): %0.2f p.p., or %0.2f' % ( 100*(models[2].eu_share[49]-models[2].eu_share[0]),
                                                                      100*(models[2].eu_share[49]/models[2].eu_share[0]-1.0) )

