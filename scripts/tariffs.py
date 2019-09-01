##########################################################################################
# This script calculates MFN tariff rates for bilateral UK-EU trade at the goods-sector level
# Results are averages of EU MFN tariffs for 6-digit HS sectors, weighted by UK-EU trade flows
# for those sectors... UK --> EU tariff is weighted by UK exports to EU; EU --> UK tariff
# is weighted by UK imports from EU. All data are for 2015.
##########################################################################################

import numpy as np
import pandas as pd

eu_names = pd.read_csv('../data/comtrade-eu-countries.csv')
eu = eu_names.country.values

wto = pd.read_csv('../data/wto-eu.csv',dtype={'hs_level':int,'hs_code':str,'avg_duty':float})
wto = wto[wto.hs_level==6]
wto.avg_duty[wto.avg_duty.isnull()] = 0.0
wto.hs_code = 'H4-'+wto.hs_code
wto=wto.drop('hs_level',axis=1)
wto=wto.rename(columns={'hs_code':'Commodity Code'})

comtrade = pd.read_csv('../data/comtrade-uk.csv')
comtrade['eu'] = comtrade.Partner.isin(eu)
comtrade = comtrade[comtrade.eu==True]
comtrade = pd.merge(left=comtrade,right=wto,how='left',on='Commodity Code')

exports = comtrade[comtrade['Trade Flow']=='Export']
imports = comtrade[comtrade['Trade Flow']=='Import']

def wavg(x,w):
    return (x*w).sum()/w.sum()


uk_to_eu_avg = wavg(exports.avg_duty,exports["Trade Value"])
eu_to_uk_avg = wavg(imports.avg_duty,imports["Trade Value"])
print '\n'
print 'Average EU MFN tariffs for 6-digit HS industries weighted by EU-UK trade flows'
print '\tUK --> EU: ' + str(uk_to_eu_avg)
print '\tEU --> UK: ' + str(eu_to_uk_avg)

print 'Scaled by goods share of total trade flows'
print '\tUK --> EU: ' + str(wavg(0.64482*exports.avg_duty,exports["Trade Value"]))
print '\tEU --> UK: ' + str(wavg(0.845511*imports.avg_duty,imports["Trade Value"]))


# percentiles

def weighted_quantile(values, quantiles, sample_weight=None, values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of initial array
    :param old_style: if True, will correct output to be consistent with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with np.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)


qs = [0.25,0.5,0.75]
q_uk_to_eu = weighted_quantile(exports.avg_duty,qs)
q_eu_to_uk = weighted_quantile(imports.avg_duty,qs)

aq_uk_to_eu = [0,0,0]
aq_uk_to_eu[0] = exports.avg_duty[exports.avg_duty<=q_uk_to_eu[0]].mean()/uk_to_eu_avg
aq_uk_to_eu[1] = exports.avg_duty[np.logical_and(exports.avg_duty>q_uk_to_eu[0],
                                                 exports.avg_duty<q_uk_to_eu[2])].mean()/uk_to_eu_avg
aq_uk_to_eu[2] = exports.avg_duty[exports.avg_duty>=q_uk_to_eu[2]].mean()/uk_to_eu_avg

aq_eu_to_uk = [0,0,0]
aq_eu_to_uk[0] = exports.avg_duty[exports.avg_duty<=q_eu_to_uk[0]].mean()/eu_to_uk_avg
aq_eu_to_uk[1] = exports.avg_duty[np.logical_and(exports.avg_duty>q_eu_to_uk[0],
                                                 exports.avg_duty<q_eu_to_uk[2])].mean()/eu_to_uk_avg
aq_eu_to_uk[2] = exports.avg_duty[exports.avg_duty>=q_eu_to_uk[2]].mean()/eu_to_uk_avg


print '\n'
print 'Average inter-quartile tariffs (relative to overall average):'
print '\t          \t0.25\t0.50\t0.75'
print '\tUK --> EU:\t%0.2f\t%0.2f\t%0.2f' % (aq_uk_to_eu[0],aq_uk_to_eu[1],aq_uk_to_eu[2])
print '\tEU --> UK:\t%0.2f\t%0.2f\t%0.2f' % (aq_eu_to_uk[0],aq_eu_to_uk[1],aq_eu_to_uk[2])
