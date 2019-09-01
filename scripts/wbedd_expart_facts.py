import pandas as pd
import numpy as np

# load the data... we need to load the CY_all dataset to know how many exporters there are total
tot = pd.read_stata('../data/CY_all.dta')
df = pd.read_stata('../data/CYD_all.dta')


tot=tot[tot.y.isin(range(2006,2009))]
df=df[df.y.isin(range(2006,2009))]

# rename columns and transform a few
tot.rename(inplace=True,columns={'A1':'nex_tot','B4i':'nd_avg'})
df.rename(inplace=True,columns={'A1':'nex',
                                'A6i':'avg_ex',
                                'B2ii':'top5_share',
                                'C2':'exit_rate',
                                'A7i':'entrant_avg_ex',
                                'A11i':'incumbent_growth_rate',
                                'A12i':'entrant_growth_rate'})

df=df[df.nex>=20]

df['entrant_exit_rate'] = 1.0-df.C3
df['entrant_rel_size'] = df.entrant_avg_ex/df.avg_ex
df['entrant_rel_growth'] = df.entrant_growth_rate - df.incumbent_growth_rate

# drop unneeded columns
tot=tot[['c','y','nex_tot','nd_avg']]
keep_cols=['c','d','y',
           'nex','avg_ex','top5_share',
           'exit_rate','entrant_exit_rate','entrant_rel_size','entrant_rel_growth']
df=df[keep_cols]

# keep only relevant countries
keep_countries=['BEL','BGR','EST','NOR','PRT','ESP','SWE']
#tot=tot[tot.c.isin(keep_countries)]
#df=df[df.c.isin(keep_countries)]

# destination grouping
eu_countries=['AUT',
              'BEL',
              'BGR',
              'CYP',
              'CZE',
              'DEN',
              'ESP',
              'EST',
              'FIN',
              'FRA',
              'GRC',
              'DEU',
              'HUN',
              'IRL',
              'ITA',
              'LVA',
              'LTU',
              'LUX',
              'MLT',
              'NLD',
              'PLD',
              'PRT',
              'ROM',
              'SVK',
              'SVN',
              'SWE']


def which_grp(d):
    if d=='GBR':
        return 'UK'
    elif d in(eu_countries):
        return 'EU'
    else:
        return 'RW'

df['cgrp'] = df.c.apply(which_grp)
df['dgrp'] = df.d.apply(which_grp)

dgrp_agg=df.groupby(['c','cgrp','dgrp','y'])['nex'].agg({'nex_max':lambda x: x.max(),
                                                         'nex_sum':lambda x: x.sum()}).reset_index()
tmp1=pd.merge(left=dgrp_agg,right=tot,how='left',on=['c','y'])
tmp1['nex_rel']=tmp1.nex_sum/tmp1.nex_tot
tmp1.loc[tmp1.dgrp!='UK','nex_rel']=tmp1.nex_rel/tmp1.nd_avg

tmp2=tmp1.drop('y',axis=1).groupby(['c','cgrp','dgrp'])['nex_rel'].mean().reset_index()
tmp3=tmp2.groupby(['cgrp','dgrp'])['nex_rel'].mean().reset_index()

tmp4=df[['top5_share','exit_rate','entrant_rel_size','entrant_exit_rate','entrant_rel_growth']].mean()
tmp5=df[np.logical_and(df.cgrp.isin(['EU','RW']),df.dgrp.isin(['UK','RW','EU']))].groupby(['dgrp','cgrp'])[['top5_share','exit_rate','entrant_rel_size','entrant_exit_rate','entrant_rel_growth']].mean()

print '\nKey facts computed from World Bank Exporter Dynamics database'

print '\nBilateral export participation rates (relative to overall export participation rate):'
print 'EU-->UK: %0.6f' % tmp3.nex_rel[np.logical_and(tmp3.cgrp=='EU',tmp3.dgrp=='UK')].values[0]
print 'EU-->RW: %0.6f' % tmp3.nex_rel[np.logical_and(tmp3.cgrp=='EU',tmp3.dgrp=='RW')].values[0]
print 'RW-->UK: %0.6f' % tmp3.nex_rel[np.logical_and(tmp3.cgrp=='RW',tmp3.dgrp=='UK')].values[0]
print 'RW-->EU: %0.6f' % tmp3.nex_rel[np.logical_and(tmp3.cgrp=='RW',tmp3.dgrp=='EU')].values[0]

print '\nDestination-level exporter distribution and dynamics:'
print 'Top 5 share\tExit rate\tEntrant rel. size\tEntrant exit rate\tEntrant rel. growth'
print '%0.6f\t%0.6f\t%0.6f\t\t%0.6f\t\t%0.6f' % (tmp4.top5_share,tmp4.exit_rate,tmp4.entrant_rel_size,tmp4.entrant_exit_rate,tmp4.entrant_rel_growth)


