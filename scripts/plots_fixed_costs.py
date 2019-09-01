#####################################################################################
#imports, small functions, etc.

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

suff = 'fixed_costs'

cnt = 0

def norm(x):
    return 100*x/x.iloc[0]

#################################################################################
# load the model results

models_uk = []
models_uk.append(pd.read_csv('../quanal/fixed_costs/output/vars_nobrexit_uk_baseline.csv'))
models_uk.append(pd.read_csv('../quanal/fixed_costs/output/vars_stoch_soft_uk_baseline.csv'))
models_uk.append(pd.read_csv('../quanal/fixed_costs/output/vars_stoch_hard_uk_baseline.csv'))
models_uk.append(pd.read_csv('../quanal/fixed_costs/output/vars_det_soft_uk_baseline.csv'))
models_uk.append(pd.read_csv('../quanal/fixed_costs/output/vars_det_hard_uk_baseline.csv'))

models_eu = []
models_eu.append(pd.read_csv('../quanal/fixed_costs/output/vars_nobrexit_eu_baseline.csv'))
models_eu.append(pd.read_csv('../quanal/fixed_costs/output/vars_stoch_soft_eu_baseline.csv'))
models_eu.append(pd.read_csv('../quanal/fixed_costs/output/vars_stoch_hard_eu_baseline.csv'))
models_eu.append(pd.read_csv('../quanal/fixed_costs/output/vars_det_soft_eu_baseline.csv'))
models_eu.append(pd.read_csv('../quanal/fixed_costs/output/vars_det_hard_eu_baseline.csv'))

models_rw = []
models_rw.append(pd.read_csv('../quanal/fixed_costs/output/vars_nobrexit_rw_baseline.csv'))
models_rw.append(pd.read_csv('../quanal/fixed_costs/output/vars_stoch_soft_rw_baseline.csv'))
models_rw.append(pd.read_csv('../quanal/fixed_costs/output/vars_stoch_hard_rw_baseline.csv'))
models_rw.append(pd.read_csv('../quanal/fixed_costs/output/vars_det_soft_rw_baseline.csv'))
models_rw.append(pd.read_csv('../quanal/fixed_costs/output/vars_det_hard_rw_baseline.csv'))

mperiods = models_uk[0].period
mseries = []

#################################################################################
# some macroeconomic plots

TREF=4
TVOTE=5
TBREXIT=8
T=20
TLR = 24
TK = 5
W=1
xlabs = [str(2011+p) for p in mperiods]
ROT=0

def autolabel(ax,rects,sign):
    for rect in rects:
        height=sign*rect.get_height()
        ax.annotate('%0.2f' % float(height),
                    (rect.get_x()+rect.get_width()/2.0,height),
                    ha='center',
                    va='center',
                    xytext=(0,sign*20),
                    textcoords='offset points')

colors=['#4daf4a','#e41a1c','#377eb8','#1b9e77','#d95f02']

def plot(all=0,expart=0):
    fig=plt.figure()
    ax=fig.add_axes([0,0,1,1])       
    
    ax.plot(mperiods[0:40],mseries[0][0:40],linestyle='-',color='black',marker=None,alpha=0.5,zorder=0,linewidth=1)
    
    if all:
        ln3=ax.plot(mperiods[(TREF-1):T],mseries[3][(TREF-1):T],dashes=(3,3),color=colors[3],linewidth=3,alpha=0.5,label='Soft (perf. foresight)')
        ln4=ax.plot(mperiods[(TREF-1):T],mseries[4][(TREF-1):T],dashes=(3,3),color=colors[4],linewidth=3,alpha=0.5,label='Hard (perf. foresight)')

    ln1=ax.plot(mperiods[(TBREXIT-1):T+1],mseries[1][(TBREXIT-1):T+1],dashes=(9,6),color=colors[0],alpha=1,label='Soft')
    ln2=ax.plot(mperiods[(TBREXIT-1):T+1],mseries[2][(TBREXIT-1):T+1],dashes=(9,6),color=colors[1],alpha=1,label='Hard')

    if expart==1:
        ln5=ax.plot(mperiods[(TREF-1):(TBREXIT+1)],mseries[2][(TREF-1):(TBREXIT+1)],linestyle='-',color=colors[2],alpha=1,marker=None,markeredgewidth=0,label='Pre-Brexit')
    else:
        ln5=ax.plot(mperiods[(TREF-1):(TBREXIT)],mseries[2][(TREF-1):(TBREXIT)],linestyle='-',color=colors[2],alpha=1,marker=None,markeredgewidth=0,label='Pre-Brexit')
                

    lns=[]
    if all:
        lns=ln1+ln2+ln5+ln3+ln4
    else:
        lns=ln1+ln2+ln5

    rects1 = ax.bar([TLR-W],[mseries[1][49]],color=colors[0],width=W,linewidth=0)
    rects2 = ax.bar([TLR+W],[mseries[2][49]],color=colors[1],width=W,linewidth=0)

    ax.set_xlim(max(TREF-TK,0),TLR+4*W)
    ax.set_xticks(range(max(TREF-TK,0),T,1),minor=True)
    ax.set_xticks(range(TREF,T,TK)+[TLR],minor=False)
    ax.set_xticklabels([xlabs[x] for x in range(TREF,T,TK)]+['Long run'],minor=False,rotation=ROT)
    sign1 = 1 if mseries[1][50]>0 else -1
    sign2 = 1 if mseries[2][50]>0 else -1

    autolabel(ax,rects1,sign1)
    autolabel(ax,rects2,sign2)

    return ax,lns

##########################
# aggregate plots

# gdp
for m in models_uk:
    mseries.append( 100*(m.rgdp/m.rgdp[0]-1.0) )

for all in [0,1]:
    plot(all)
    plt.ylim(-2,0.5)
    plt.savefig('fig/gdp_'+suff+str(all)+'.pdf')
    plt.clf()

# investment
mseries = []
for m in models_uk:
    mseries.append( 100*((m.x/m.ngdp)/(m.x[0]/m.ngdp[0])-1.0) )
    
for all in [0,1]:
    plot(all)
    plt.ylim(-3,2.0)
    #plt.yticks([-4,-3,-2,-1,0,1],['-4.0','-3.0','-2.0','-1.0','0.0','1.0'])
    plt.savefig('fig/inv_'+suff+str(all)+'.pdf')
    plt.clf()
    plt.close()

# consumption
mseries = []
for m in models_uk:
    mseries.append( 100*((m.c)/(m.c[0])-1.0) )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-2,0.5)
    plt.legend(lns,labs,prop={'size':20},loc='lower left',ncol=2)
    plt.savefig('fig/cons_'+suff+str(all)+'.pdf')
    plt.clf()
    plt.close()
    
# im
mseries = []
for m in models_uk:
    #mseries.append( 100*((m.rim1+m.rim2)/(m.rim1[0]+m.rim2[0])-1.0) )
    mseries.append( 100*( ((m.rim1+m.rim2)/m.ngdp)/ ((m.rim1[0]+m.rim2[0])/m.ngdp[0]) -1.0 ) )

for all in [0,1]:
    plot(all)
    #plt.ylim(-25,5)
    plt.savefig('fig/im_'+suff+str(all)+'.pdf')
    plt.clf()
    plt.close()

##########################
# trade with EU plots

# exports
mseries = []
for m in models_uk:
    mseries.append( 100*((m.rex1/m.ngdp)/(m.rex1[0]/m.ngdp[0])-1.0) )
    
for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-60,10)
    #plt.yticks([-60,-50,-40,-30,-20,-10,0,10],['-60.0','-50.0','-40.0','-30.0','-20.0','-10.0','0.0','10.0'])
    plt.legend(lns,labs,prop={'size':20},loc='lower left',ncol=2)
    plt.savefig('fig/eu_exports_'+suff+str(all)+'.pdf')
    plt.clf()

# imports
mseries = []
for m in models_uk:
    mseries.append( 100*((m.rim1/m.ngdp)/(m.rim1[0]/m.ngdp[0])-1.0) )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-60,10)
    #plt.yticks([-60,-50,-40,-30,-20,-10,0,10],['-60.0','-50.0','-40.0','-30.0','-20.0','-10.0','0.0','10.0'])
    plt.savefig('fig/eu_imports_'+suff+str(all)+'.pdf')
    plt.clf()

# export participation
mseries = []
for m in models_uk:
    mseries.append( 100*(m.exrate1/m.exrate1[0]-1.0) )

for all in [0,1]:
    ax,lns=plot(all,1)
    labs=[l.get_label() for l in lns]
    plt.ylim(-60,10)
    plt.savefig('fig/eu_expart_'+suff+str(all)+'.pdf')
    plt.clf()
    
# import participation
mseries = []
cnt=0
for m in models_eu:
    mseries.append( 100*(m.exrate0/m.exrate0[0]-1.0) )
    cnt = cnt+1

for all in [0,1]:
    ax,lns=plot(all,1)
    labs=[l.get_label() for l in lns]
    plt.ylim(-60,10)
    plt.savefig('fig/eu_impart_'+suff+str(all)+'.pdf')
    plt.clf()

plt.close('all')
