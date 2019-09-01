##########################################################################################
# This script calculates MFN tariff rates for bilateral UK-EU trade at the goods-sector level
# Results are averages of EU MFN tariffs for 6-digit HS sectors, weighted by UK-EU trade flows
# for those sectors... UK --> EU tariff is weighted by UK exports to EU; EU --> UK tariff
# is weighted by UK imports from EU. All data are for 2015.
##########################################################################################

import numpy as np
import pandas as pd

ntb = pd.read_csv('../../data/ntb.csv')

data = pd.read_pickle('wiod_uk_ts.pik')
data = data[data.year==2011]
data = data[data.region=='EU']
data = pd.merge(left=data,right=ntb,how='left',on='industry_code')
#goods = data[data.sector=='goods']
#svcs = data[data.sector=='services']

def wavg(x,w):
    return (x*w).sum()/w.sum()

print '\n'
print 'Average EU-USA reducible NTBs for WIOD industries weighted by EU-UK trade flows'
print 'UK imports from EU:' + str(wavg(data.ntb,data.total_imports))
#print 'Intermediate goods:    ' + str(wavg(goods.ntb,goods.intermediate_imports))
#print 'Final goods:           ' + str(wavg(goods.ntb,goods.final_imports))
#print 'Intermediate services: ' + str(wavg(svcs.ntb,svcs.intermediate_imports))
#print 'Final services:        ' + str(wavg(svcs.ntb,svcs.final_imports))
print 'EU imports from UK:' + str(wavg(data.ntb,data.total_exports))
#print 'Intermediate goods:    ' + str(wavg(goods.ntb,goods.intermediate_exports))
#print 'Final goods:           ' + str(wavg(goods.ntb,goods.final_exports))
#print 'Intermediate services: ' + str(wavg(svcs.ntb,svcs.intermediate_exports))
#print 'Final services:        ' + str(wavg(svcs.ntb,svcs.final_exports))
