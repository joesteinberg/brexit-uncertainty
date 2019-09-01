import pandas as pd
import numpy as np

df = pd.read_csv('../data/EFIGE_Truncated_1_updated.csv',sep=',')
df=df[df.country.isin(['FRA','GER','ITA','SPA','UK'])]

df['d4'] = df.d4.fillna(0.0)
df['d13_1']=df.d13_1.fillna(0.0)
df['d13_2']=df.d13_2.fillna(0.0)
df['d13_3']=df.d13_3.fillna(0.0)

##########################################################################################
# overall export participation rates

df['exporter']=df.d4>1.0e-6
expart_rates = df.groupby('country').agg({'mark':lambda x: x.nunique(),'exporter':np.mean}).reset_index()

##########################################################################################
# how many export to the EU, rest of the world

df['eu_export_frac']=df.d13_1+df.d13_2
df['rw_export_frac']=100.0-df.eu_export_frac
df['eu']=df.eu_export_frac>1.0e-6
df['rw']=df.rw_export_frac>1.0e-6

expart_rates_by_dest = df[df.exporter].groupby('country')[['eu','rw']].mean().reset_index()



print '\nKey facts computed from EFIGE database...'

print'\nOverall export participation_rates:'
print 'UK: %0.6f' % expart_rates.exporter[expart_rates.country=='UK'].values[0]
print 'EU: %0.6f' % expart_rates.exporter[expart_rates.country!='UK'].mean()

print'\nBilateral export participation_rates (conditional on exporting):'
print 'UK-->EU: %0.6f' % expart_rates_by_dest.eu[expart_rates_by_dest.country=='UK'].values[0]
print 'UK-->RW: %0.6f' % expart_rates_by_dest.rw[expart_rates_by_dest.country=='UK'].values[0]
print 'EU-->RW: %0.6f' % expart_rates_by_dest.rw[expart_rates_by_dest.country!='UK'].mean()
