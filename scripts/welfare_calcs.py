import numpy as np
import pandas as pd
import locale
locale.setlocale(locale.LC_ALL,'en_US.utf8')

path = '../quanal/dyn_mkt_pen/output/'
path_ac = '../quanal/fixed_costs/output/'
path_multi = '../quanal/multisector/output/'

######################################################################################################
# Welfare analysis for Brexit exercise

T2015 = 4
T2016 = 5
TBGP = 100

#phi = 0.3957
phi=1.0

S=['baseline',
   'dynamic_sunk_cost',
   'static_arkolakis',
   'static_melitz',
   'no_export_costs',
   'ac',
   'highpi',
   'lowpi',
   'kappa_uncertainty',
   'kappa_ntb_uncertainty',
   'rb_perm',
   'rb_temp',
   'idio_uncertainty',
   'long_prebrexit',
   'multi',
   'multif',
   'finaut',
   'szeta',
   'spsi',
   'sxr',
   'sld',
   'shighr']

N=['Baseline',
   'Dynamic sunk cost',
   'Static market. pen.',
   'Static fixed cost',
   'No export costs',
   'Alessandria-Choi',
   'Lower prob. of hard Brexit',
   'Higher prob. of hard Brexit',
   'Marketing costs, not icebergs',
   'Marketing costs and icebergs',
   'Reversible Brexit (permanent)',
   'Reversible Brexit (temporary)',
   'Firm-level policy uncertainty',
   'Longer pre-Brexit period',
   'Multi-sector',
   'Multi-sector w/ frictions',
   'Financial autarky',
   'Lower Armington elasticity',
   'Higher risk aversion',
   'Lower exit rate',
   'Lower customer depreciation',
   'Higher interest rate']

def welfare(suff):
    models=[]
    if(suff=='multi'):
        models.append(pd.read_csv(path_multi+'vars_nobrexit_uk.csv'))
        models.append(pd.read_csv(path_multi+'vars_stoch_opt_uk.csv'))
        models.append(pd.read_csv(path_multi+'vars_stoch_pes_uk.csv'))
        models.append(pd.read_csv(path_multi+'vars_det_opt_uk.csv'))
        models.append(pd.read_csv(path_multi+'vars_det_pes_uk.csv'))
    elif(suff=='multif'):
        models.append(pd.read_csv(path_multi+'vars_nobrexit_uk_frictions.csv'))
        models.append(pd.read_csv(path_multi+'vars_stoch_opt_uk_frictions.csv'))
        models.append(pd.read_csv(path_multi+'vars_stoch_pes_uk_frictions.csv'))
        models.append(pd.read_csv(path_multi+'vars_det_opt_uk_frictions.csv'))
        models.append(pd.read_csv(path_multi+'vars_det_pes_uk_frictions.csv'))
    elif(suff=='ac'):
        models.append(pd.read_csv(path_ac+'vars_nobrexit_uk_baseline.csv'))
        models.append(pd.read_csv(path_ac+'vars_stoch_soft_uk_baseline.csv'))
        models.append(pd.read_csv(path_ac+'vars_stoch_hard_uk_baseline.csv'))
        models.append(pd.read_csv(path_ac+'vars_det_soft_uk_baseline.csv'))
        models.append(pd.read_csv(path_ac+'vars_det_hard_uk_baseline.csv'))
    elif(suff=='rb_perm'):
        models.append(pd.read_csv(path+'vars_nobrexit_uk_baseline.csv'))
        models.append(pd.read_csv(path+'vars_stoch_rb_perm_soft_uk_baseline.csv'))
        models.append(pd.read_csv(path+'vars_stoch_rb_perm_hard_uk_baseline.csv'))
        models.append(pd.read_csv(path+'vars_det_soft_uk_baseline.csv'))
        models.append(pd.read_csv(path+'vars_det_hard_uk_baseline.csv'))
    elif(suff=='rb_temp'):
        models.append(pd.read_csv(path+'vars_nobrexit_uk_baseline.csv'))
        models.append(pd.read_csv(path+'vars_stoch_rb_temp_soft_uk_baseline.csv'))
        models.append(pd.read_csv(path+'vars_stoch_rb_temp_hard_uk_baseline.csv'))
        models.append(pd.read_csv(path+'vars_det_temp_soft_uk_baseline.csv'))
        models.append(pd.read_csv(path+'vars_det_temp_hard_uk_baseline.csv'))
    else:
        models.append(pd.read_csv(path+'vars_nobrexit_uk_'+suff+'.csv'))
        models.append(pd.read_csv(path+'vars_stoch_soft_uk_'+suff+'.csv'))
        models.append(pd.read_csv(path+'vars_stoch_hard_uk_'+suff+'.csv'))
        models.append(pd.read_csv(path+'vars_det_soft_uk_'+suff+'.csv'))
        models.append(pd.read_csv(path+'vars_det_hard_uk_'+suff+'.csv'))

    w0_2015 = models[0].W[T2015]
    wct = []
    wct.append((w0_2015/models[1].W[T2015])**(1/phi)-1.0)
    wct.append((w0_2015/models[2].W[T2015])**(1/phi)-1.0)

    wclr = []
    wclr.append((models[0].c[0]/models[1].c.tail(1).values[0]-1.0))
    wclr.append((models[0].c[0]/models[2].c.tail(1).values[0]-1.0))

    w3_2015 = models[3].W[T2015]
    w4_2015 = models[4].W[T2015]
    wcu = []
    wcu.append((w3_2015/models[1].W[T2015])**(1/phi)-1.0)
    wcu.append((w4_2015/models[2].W[T2015])**(1/phi)-1.0)

    wcd=[]
    for i in range(len(wcu)):
        wcu[i] = abs(100*wcu[i]/wct[i])
        wcd.append(100*(1.0-wclr[i]/wct[i]))
        wct[i] = 100*wct[i]  
        wclr[i]=wclr[i]*100
    
    strs=[locale.format('%0.2f',x,grouping=True) for x in wct+wcu]

    return strs


# open the tex file
file=open('tex/results-welfare.tex','wb')

file.write('\\begin{table}[h!]\n')
file.write('\\renewcommand{\\arraystretch}{1.2}\n')
file.write('\\begin{center}\n')
file.write('\\caption{U.K. welfare losses from Brexit}\n')
file.write('\\label{tab:results-welfare}\n')

# number of columns
file.write('\\begin{tabular}{lcccc}')
file.write('\\toprule\n')

# headers
file.write('&\\multicolumn{2}{c}{Total (cons. equiv.)}&\\multicolumn{2}{c}{Uncertainty (pct. total)}\\\\\n')
file.write('\\cmidrule(rl){2-3}\\cmidrule(rl){4-5}\n')
file.write('Model')
file.write('&\\multicolumn{1}{p{1.5cm}}{\\centering Soft}&\\multicolumn{1}{p{1.5cm}}{\\centering Hard}')
file.write('&\\multicolumn{1}{p{1.5cm}}{\\centering Soft}&\\multicolumn{1}{p{1.5cm}}{\\centering Hard}')
file.write('\\\\\n')
file.write('\\midrule\n')

# baseline model
suff=S[0]
strs=welfare(suff)
file.write(N[0])
for s in strs:
    file.write('&'+s)  
file.write('\\\\\n')
file.write('\\\\\n')

# alternative scenarios
file.write('\\multicolumn{5}{l}{\\textit{(a) Alternative scenarios}}\\\\\n')
for i in range(6,14):
    suff=S[i]
    strs=welfare(suff)
    file.write(N[i])
    for s in strs:
        file.write('&'+s)  
    file.write('\\\\\n')
file.write('\\\\\n')

# alternative model setup
file.write('\\multicolumn{5}{l}{\\textit{(b) Alternative models}}\\\\\n')

for i in range(1,6):
    suff=S[i]
    strs=welfare(suff)
    file.write(N[i])
    for s in strs:
        file.write('&'+s)  
    file.write('\\\\\n')
file.write('\\\\\n')

# sensitivity analyses
file.write('\\multicolumn{5}{l}{\\textit{(c) Sensitivity analyses}}\\\\\n')
for i in range(14,22):
    suff=S[i]
    strs=welfare(suff)
    file.write(N[i])
    for s in strs:
        file.write('&'+s)  
    file.write('\\\\\n')


file.write('\\bottomrule\n')
file.write('\\end{tabular}\n')
file.write('\\end{center}\n')
file.write('\\end{table}\n')
#file.write('\\end{landscape}\n')
file.close()
