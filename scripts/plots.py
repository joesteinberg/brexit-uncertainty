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

suff = 'baseline'

cnt = 0

def norm(x):
    return 100*x/x.iloc[0]

#################################################################################
# load the model results

models_uk = []
models_uk.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_nobrexit_uk_'+suff+'.csv'))
models_uk.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_stoch_soft_uk_'+suff+'.csv'))
models_uk.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_stoch_hard_uk_'+suff+'.csv'))
models_uk.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_det_soft_uk_'+suff+'.csv'))
models_uk.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_det_hard_uk_'+suff+'.csv'))

models_eu = []
models_eu.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_nobrexit_eu_'+suff+'.csv'))
models_eu.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_stoch_soft_eu_'+suff+'.csv'))
models_eu.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_stoch_hard_eu_'+suff+'.csv'))
models_eu.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_det_soft_eu_'+suff+'.csv'))
models_eu.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_det_hard_eu_'+suff+'.csv'))

models_rw = []
models_rw.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_nobrexit_rw_'+suff+'.csv'))
models_rw.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_stoch_soft_rw_'+suff+'.csv'))
models_rw.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_stoch_hard_rw_'+suff+'.csv'))
models_rw.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_det_soft_rw_'+suff+'.csv'))
models_rw.append(pd.read_csv('../quanal/dyn_mkt_pen/output/vars_det_hard_rw_'+suff+'.csv'))

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

models_uk[1].texrate1[0:TBREXIT] = models_uk[1].exrate1[0:TBREXIT]

models_uk[1].texrate2[0:TBREXIT] = 0.5*(models_uk[3].texrate2[0:TBREXIT] + models_uk[4].texrate2[0:TBREXIT])
models_rw[1].texrate0[0:TBREXIT] = 0.5*(models_rw[3].texrate0[0:TBREXIT] + models_rw[4].texrate0[0:TBREXIT])
models_uk[2].texrate2[0:TBREXIT] = 0.5*(models_uk[3].texrate2[0:TBREXIT] + models_uk[4].texrate2[0:TBREXIT])
models_rw[2].texrate0[0:TBREXIT] = 0.5*(models_rw[3].texrate0[0:TBREXIT] + models_rw[4].texrate0[0:TBREXIT])

models_uk[1].mktpen2[0:TBREXIT] = 0.5*(models_uk[3].mktpen2[0:TBREXIT] + models_uk[4].mktpen2[0:TBREXIT])
models_rw[1].mktpen0[0:TBREXIT] = 0.5*(models_rw[3].mktpen0[0:TBREXIT] + models_rw[4].mktpen0[0:TBREXIT])
models_uk[2].mktpen2[0:TBREXIT] = 0.5*(models_uk[3].mktpen2[0:TBREXIT] + models_uk[4].mktpen2[0:TBREXIT])

#models_rw[2].mktpen0[

models_rw[2].mktpen0[0:TBREXIT] = 0.5*(models_rw[3].mktpen0[0:TBREXIT] + models_rw[4].mktpen0[0:TBREXIT])



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

def plot(all=0):
    fig=plt.figure()
    ax=fig.add_axes([0,0,1,1])       
    
    ax.plot(mperiods[0:40],mseries[0][0:40],linestyle='-',color='black',marker=None,alpha=0.5,zorder=0,linewidth=1)
    
    if all:
        ln3=ax.plot(mperiods[(TREF-1):T],mseries[3][(TREF-1):T],dashes=(3,3),color=colors[3],linewidth=3,alpha=0.5,label='Soft (perf. foresight)')
        ln4=ax.plot(mperiods[(TREF-1):T],mseries[4][(TREF-1):T],dashes=(3,3),color=colors[4],linewidth=3,alpha=0.5,label='Hard (perf. foresight)')

    ln1=ax.plot(mperiods[(TBREXIT-1):T],mseries[1][(TBREXIT-1):T],dashes=(9,6),color=colors[0],alpha=1,label='Soft')
    ln2=ax.plot(mperiods[(TBREXIT-1):T],mseries[2][(TBREXIT-1):T],dashes=(9,6),color=colors[1],alpha=1,label='Hard')

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
    if(all==1):
        plt.ylabel('pct. change')
        plt.title('Impact of Brexit on UK consumption',y=1.05)
        plt.savefig('fig/cons_web.png')
    plt.clf()
    plt.close()

# ex
mseries = []
for m in models_uk:
    mseries.append( 100*( ((m.ex1+m.ex2)/m.ngdp)/ ((m.ex1[0]+m.ex2[0])/m.ngdp[0]) -1.0 ) )

for all in [0,1]:
    plot(all)
    #plt.ylim(-25,5)
    plt.savefig('fig/ex_'+suff+str(all)+'.pdf')
    plt.clf()
    plt.close()

# im
mseries = []
for m in models_uk:
    #mseries.append( 100*((m.rim1+m.rim2)/(m.rim1[0]+m.rim2[0])-1.0) )
    mseries.append( 100*( ((m.im1+m.im2)/m.ngdp)/ ((m.im1[0]+m.im2[0])/m.ngdp[0]) -1.0 ) )

for all in [0,1]:
    plot(all)
    #plt.ylim(-25,5)
    plt.savefig('fig/im_'+suff+str(all)+'.pdf')
    plt.clf()
    plt.close()

# nx
mseries = []
for m in models_uk:
    mseries.append( 100*(m.nx1+m.nx2)/m.ngdp - 100*(m.nx1[0]+m.nx2[0])/m.ngdp[0] )

for all in [0,1]:
    plot(all)
    plt.ylim(-.4,0.8)
    #plt.yticks([-.2,0.0,0.2,.4,.6,.8],['-0.2','0.0','0.2','0.4','0.6','0.8'])
    plt.savefig('fig/nx_'+suff+str(all)+'.pdf')
    plt.clf()
    plt.close()

##########################
# trade with EU plots

# exports
mseries = []
for m in models_uk:
    mseries.append( 100*((m.ex1/m.ngdp)/(m.ex1[0]/m.ngdp[0])-1.0) )
    
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
    mseries.append( 100*((m.im1/m.ngdp)/(m.im1[0]/m.ngdp[0])-1.0) )
    #mseries.append( 100*(m.rim1/m.rim1[0]-1.0) )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-60,10)
    #plt.yticks([-60,-50,-40,-30,-20,-10,0,10],['-60.0','-50.0','-40.0','-30.0','-20.0','-10.0','0.0','10.0'])
    plt.savefig('fig/eu_imports_'+suff+str(all)+'.pdf')
    plt.clf()

# export participation
mseries = []
cnt=0
for m in models_uk:
    if cnt==1 or cnt==3:
        mseries.append( 100*(m.texrate1/m.texrate1[0]-1.0) )
    else:
        mseries.append( 100*(m.exrate1/m.exrate1[0]-1.0) )
    cnt = cnt+1

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-50,10)
    plt.savefig('fig/eu_expart_'+suff+str(all)+'.pdf')
    plt.clf()

# export mkt pen
mseries = []
for m in models_uk:
    mseries.append( 100*(m.mktpen1/m.mktpen1[0]-1.0) )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-30,5)
    plt.savefig('fig/eu_ex_mktpen_'+suff+str(all)+'.pdf')
    plt.clf()
    
# import participation
mseries = []
cnt=0
for m in models_eu:
    #if cnt==1 or cnt==3:
    #    mseries.append( 100*(m.texrate0/m.texrate0[0]-1.0) )
    #else:
    mseries.append( 100*(m.exrate0/m.exrate0[0]-1.0) )
    cnt = cnt+1

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-50,10)
    plt.savefig('fig/eu_impart_'+suff+str(all)+'.pdf')
    plt.clf()

# import mkt pen
mseries = []
for m in models_eu:
    mseries.append( 100*(m.mktpen0/m.mktpen0[0]-1.0) )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-30,5)
    plt.savefig('fig/eu_im_mktpen_'+suff+str(all)+'.pdf')
    plt.clf()

# net exports
mseries = []
for m in models_uk:
    mseries.append( 100*m.nx1/m.ngdp - 100*m.nx1[0]/m.ngdp[0] )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-0.25,2.0)
    #plt.yticks([-0.2,0.0,0.2,0.4,0.6,0.8,1.0],['-0.2','0.0','0.2','0.4','0.6','0.8','1.0'])
    plt.savefig('fig/eu_nx_'+suff+str(all)+'.pdf')
    plt.clf()

# rer
mseries = []
for m in models_uk:
    mseries.append( 100*(m.rer1/m.rer1[0]-1.0) )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-1.5,0.25)
    #plt.yticks([-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2],['-1.4','-1.2','-1.0','-0.8','-0.6','-0.4','-0.2','0.0','0.2'])
    plt.savefig('fig/eu_rer_'+suff+str(all)+'.pdf')
    plt.clf()

##########################
# trade with RW plots

# exports
mseries = []
for m in models_uk:
    mseries.append( 100*((m.ex2/m.ngdp)/(m.ex2[0]/m.ngdp[0])-1.0) )
    #mseries.append( 100*(m.rex2/m.rex2[0]-1.0) )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-3,10)
    #plt.yticks([-2,0,2,4,6,8,10,12],['-2.0','0.0','2.0','4.0','6.0','8.0','10.0','12.0'])
    plt.legend(lns,labs,prop={'size':20},loc='upper left',ncol=2)
    plt.savefig('fig/rw_exports_'+suff+str(all)+'.pdf')
    plt.clf()

# imports
mseries = []
for m in models_uk:
    mseries.append( 100*((m.im2/m.ngdp)/(m.im2[0]/m.ngdp[0])-1.0) )
    #mseries.append( 100*(m.rim2/m.rim2[0]-1.0) )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-3,10)
    #plt.yticks([-2,0,2,4,6,8,10],['-2.0','0.0','2.0','4.0','6.0','8.0','10.0','12.0'])
    plt.savefig('fig/rw_imports_'+suff+str(all)+'.pdf')
    plt.clf()

# export participation
mseries = []
for m in models_uk:
    mseries.append( 100*(m.texrate2/m.texrate2[0]-1.0) )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-2,6)
    plt.savefig('fig/rw_expart_'+suff+str(all)+'.pdf')
    plt.clf()

# export mkt pen
mseries = []
for m in models_uk:
    mseries.append( 100*(m.mktpen2/m.mktpen2[0]-1.0) )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-3,10)
    #plt.yticks([-4,0,4,8,12,16,20],['-4.0','0.0','4.0','8.0','12.0','16.0','20.0'])
    plt.savefig('fig/rw_ex_mktpen_'+suff+str(all)+'.pdf')
    plt.clf()

# export participation
mseries = []
for m in models_rw:
    mseries.append( 100*(m.texrate0/m.texrate0[0]-1.0) )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-2,6)
    plt.savefig('fig/rw_impart_'+suff+str(all)+'.pdf')
    plt.clf()

# import market penetration
mseries = []
for m in models_rw:
    mseries.append( 100*(m.mktpen0/m.mktpen0[0]-1.0) )
    
for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-3,10)
    #plt.yticks([-4,0,4,8,12,16,20],['-4.0','0.0','4.0','8.0','12.0','16.0','20.0'])
    plt.savefig('fig/rw_im_mktpen_'+suff+str(all)+'.pdf')
    plt.clf()

# net exports
mseries = []
for m in models_uk:
    mseries.append( 100*m.nx2/m.ngdp - 100*m.nx2[0]/m.ngdp[0] )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-2,0.5)
    #plt.yticks([-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4],['-1.0','-0.8','-0.6','-0.4','-0.2','0.0','0.2','0.4'])
    plt.savefig('fig/rw_nx_'+suff+str(all)+'.pdf')
    plt.clf()

# rer
mseries = []
for m in models_uk:
    mseries.append( 100*(m.rer2/m.rer2[0]-1.0) )

for all in [0,1]:
    ax,lns=plot(all)
    labs=[l.get_label() for l in lns]
    plt.ylim(-1.5,0.25)
    #plt.yticks([-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2],['-1.4','-1.2','-1.0','-0.8','-0.6','-0.4','-0.2','0.0','0.2'])
    plt.savefig('fig/rw_rer_'+suff+str(all)+'.pdf')
    plt.clf()
    
    plt.close('all')
