##########################################################################################
# This script loads raw WIOD data and creates the input files used by MATLAB
# to calibrate the model, and by Python to create alternative IO matrices
##########################################################################################

##########################################################################################
# Imports, constants, etc.
##########################################################################################

import pandas as pd
import numpy as np
import regions

sector = 0

##########################################################################################
# Main script body
##########################################################################################

work=True

def prepare_data(home,year=2011):

    # pull all records where the home country is the destination
    filter = ""

    if home == 'Constructed RoW':
        filter = "("
        first=True
        for region in regions.regions.keys():
            for country in regions.regions[region]:
                if first==True:
                    filter = filter + "('in_country'!=%s" % country + ")"
                else:
                    filter = filter + "&('in_country'!=%s" % country + ")"
                first=False
        filter = filter + ")"

    else:
        filter = "("
        first=True
        for country in regions.regions[home]:
            if first==True:
                filter = filter + "('in_country'==%s" % country + ")"
            else:
                filter = filter + "|('in_country'==%s" % country + ")"
            first=False
        filter = filter + ")"

    filter = filter + "&('year'==%s" % str(year) + ")"

    print '\nQuery string: '+filter+'\n'

    try:
        selected = store.select('df',filter)
        print 'Number of rows: '+str(len(selected)) + '\n'
    except Exception, e:
        print 'Error selecting data from HFDStore: %s' % e + '\n'
        print 'Query: ' + filter
        return []
    else:
        # add in regions data
        selected['in_region'] = selected['in_country'].apply(lambda s: regions.which_region(s))
        selected['out_region'] = selected['out_country'].apply(lambda s: regions.which_region(s))

        # mask home as destination
        if(home != 'Constructed RoW'):
            mask_in = selected['in_region'] == home
        else:
            mask_in = selected['in_country'] != ''
            for region in regions.regions.keys():
                mask_in = np.logical_and(mask_in,selected['in_region'] != region)

        # gross output and value added .....................................................

        # mask for value added in home country
        mask_VA = np.logical_or(selected['out_industry_code']=='r99',
                                selected['out_industry_code']=='r61')
        mask_VA = np.logical_or(mask_VA,
                                selected['out_industry_code']=='r62')
        mask_VA = np.logical_or(mask_VA,
                                selected['out_industry_code']=='r63')
        mask_VA = np.logical_or(mask_VA,
                                selected['out_industry_code']=='r64')
        mask_VA = np.logical_or(mask_VA,
                                selected['out_industry_name']=='International Transport Margins')

        mask_VA = np.logical_and(mask_VA,selected['in_type']=='I')
        mask_VA = np.logical_and(mask_VA,mask_in)
                                                 
        # mask for gross output in home country (calculated as total intermediate
        # usage plus value added)
        mask_GO = np.logical_and(mask_in,selected['in_type']=='I')
        mask_GO = np.logical_and(mask_GO,np.logical_or(
                selected['out_industry_code']=='r64',
                selected['out_industry_code']=='r60'))

        # domestic gross output and value added
        #grouped = selected[mask_VA].groupby('s')
        va = selected[mask_VA]['value'].sum()

        #grouped = selected[mask_GO].groupby('s')
        #go = grouped['value'].agg(np.sum)
        va_df = pd.DataFrame({'va':[va]})
        va_df['i']=regions.regions_dict[home]

        # intermediate inputs and final uses ..............................................

        # mask for intermediate inputs into home country
        # NOTE: we want both international and intranational transactions
        # so as to capture imported and domestic intermediates
        mask_inputs = np.logical_and(mask_in,selected['in_type']=='I')
        mask_inputs = np.logical_and(mask_inputs,selected['out_type']=='I')
        mask_inputs = np.logical_and(mask_inputs,selected['out_country']!='')
        mask_inputs = np.logical_and(mask_inputs,selected['out_country']!='TOT')

        # mask for final uses in home country
        mask_CONS = np.logical_and(mask_in,selected['in_type']=='F')
        mask_CONS = np.logical_and(mask_CONS,selected['out_type']=='I')
        mask_CONS = np.logical_and(mask_CONS,selected['out_country']!='')
        mask_CONS = np.logical_and(mask_CONS,selected['out_country']!='TOT')
        mask_CONS = np.logical_and(mask_CONS,
                                   np.logical_or(selected['in_industry_code']=='c37',
                                                 np.logical_or(selected['in_industry_code']=='c38',
                                                               selected['in_industry_code']=='c39')))

        mask_INV = np.logical_and(mask_in,selected['in_type']=='F')
        mask_INV = np.logical_and(mask_INV,selected['out_type']=='I')
        mask_INV = np.logical_and(mask_INV,selected['out_country']!='')
        mask_INV = np.logical_and(mask_INV,selected['out_country']!='TOT')
        mask_INV = np.logical_and(mask_INV,
                                   np.logical_or(selected['in_industry_code']=='c41',
                                                 selected['in_industry_code']=='c42'))

        # now aggregate inputs from all countries in list
        inin_df=0
        cons_df=0
        inv_df=0
        
        first=True
        for source in regions.regions.keys():               
            mask_inputs2 = np.logical_and(mask_inputs,
                                          selected['out_region']==source)                
            mask_cons2 = np.logical_and(mask_CONS,
                                        selected['out_region']==source)
            mask_inv2 = np.logical_and(mask_INV,
                                       selected['out_region']==source)
                
            #grouped = selected[mask_inputs2].groupby(['s','r'])
            inin2 = selected[mask_inputs2]['value'].sum()
            inin2_df = pd.DataFrame({'m':[inin2]})
            inin2_df['i']=regions.regions_dict[home]
            inin2_df['j']=regions.regions_dict[source]
            
            #grouped = selected[mask_cons2].groupby('r')
            cc2 = selected[mask_cons2]['value'].sum()
            cc2_df = pd.DataFrame({'c':[cc2]})
            cc2_df['i']=regions.regions_dict[home]
            cc2_df['j']=regions.regions_dict[source]
                
            #grouped = selected[mask_inv2].groupby('r')
            inv2 = selected[mask_inv2]['value'].sum()       
            inv2_df = pd.DataFrame({'x':[inv2]})
            inv2_df['i']=regions.regions_dict[home]
            inv2_df['j']=regions.regions_dict[source]
                
            if first:
                inin_df=inin2_df
                cc_df=cc2_df
                inv_df=inv2_df
            else:
                inin_df = inin_df.append(inin2_df)
                cc_df = cc_df.append(cc2_df)
                inv_df = inv_df.append(inv2_df)
            
            first=False
                
        # now do aggregate across all countries other than those in user-supplied
        # list to get RoW
        mask_inputs2 = mask_inputs
        mask_cons2 = mask_CONS
        mask_inv2 = mask_INV
        for source in regions.regions.keys():
            mask_inputs2 = np.logical_and(mask_inputs2,selected['out_region']!=source)
            mask_cons2 = np.logical_and(mask_cons2,selected['out_region']!=source)
            mask_inv2 = np.logical_and(mask_inv2,selected['out_region']!=source)

        #grouped = selected[mask_inputs2].groupby(['s','r'])
        inin2 = selected[mask_inputs2]['value'].sum()
        inin2_df = pd.DataFrame({'m':[inin2]})
        inin2_df['i']=regions.regions_dict[home]
        inin2_df['j']=regions.regions_dict['Constructed RoW']
        #inin2_df['j']=regions.nr
        inin_df = inin_df.append(inin2_df)

        #grouped = selected[mask_cons2].groupby('r')
        cc2 = selected[mask_cons2]['value'].sum()
        cc2_df = pd.DataFrame({'c':[cc2]})
        cc2_df['i']=regions.regions_dict[home]
        cc2_df['j']=regions.regions_dict['Constructed RoW']
        #cc2_df['j']=regions.nr
        cc_df = cc_df.append(cc2_df)

        #grouped = selected[mask_inv2].groupby('r')
        inv2 = selected[mask_inv2]['value'].sum()
        inv2_df = pd.DataFrame({'x':[inv2]})
        inv2_df['i']=regions.regions_dict[home]
        inv2_df['j']=regions.regions_dict['Constructed RoW']
        #inv2_df['j']=regions.nr
        inv_df = inv_df.append(inv2_df)

        va_df=va_df.reset_index()
        inin_df=inin_df.reset_index()
        cc_df=cc_df.reset_index()
        inv_df=inv_df.reset_index()

        fin_df = pd.merge(left=cc_df,right=inv_df,how='left',on=['i','j'])
        fin_df['f'] = fin_df['c']+fin_df['x']
        fin_df = fin_df.drop(['c','x'],axis=1)

        return va_df, inin_df, fin_df

if work:
    
    # open HDFStore
    store = pd.HDFStore('/home/joe/Datasets/WIOD/hdf/wiot.h5','r')
    print 'store:\n%s' % store + '\n'

    # set up list of countries
    regions2 = regions.regions.keys() + ['Constructed RoW']

    years=[2011]
    for year in years:
        print '\nProcessing data for ' + str(year) + '...'
        va_df=0
        inin_df=0
        fin_df=0
        first=True
        for region in regions2:
            print '\t...' + region
            va2_df,inin2_df,fin2_df=prepare_data(region,year)

            if first==True:
                va_df = va2_df
                inin_df=inin2_df
                fin_df=fin2_df
            else:
                va_df = va_df.append(va2_df)
                inin_df = inin_df.append(inin2_df)
                fin_df = fin_df.append(fin2_df)
            first=False
                
        tmp = inin_df.groupby(['i'])['m'].sum()
        tmp = tmp.reset_index()
        va_df = pd.merge(left=va_df,right=tmp,how='left',on = ['i'])
        va_df['go'] = va_df['va']+va_df['m']

        va_df.to_csv(regions.path + 'va' + str(year) + '.txt',
                     sep=' ',columns=['i','go','va'],index=False)

        inin_df.to_csv(regions.path + 'inin'+str(year) + '.txt',
                       sep=' ',columns=['i','j','m'],index=False)

        fin_df.to_csv(regions.path + 'fin'+str(year) + '.txt',
                      sep=' ',columns=['i','j','f'],index=False)

        with open(regions.path + 'va' + str(year) + '.txt', 'a') as myfile:
            myfile.write('\n')

        with open(regions.path + 'inin' + str(year) + '.txt', 'a') as myfile:
            myfile.write('\n')

        with open(regions.path + 'fin' + str(year) + '.txt', 'a') as myfile:
            myfile.write('\n')

    store.close()

    # check market clearing: gross output (which equals value added plus intermediates by definition)
    # should equal intermediate use plus final use
    totalm = inin_df.groupby(['j'])['m'].sum().reset_index().rename(columns={'m':'md','j':'i'})
    totalf = fin_df.groupby(['j'])['f'].sum().reset_index().rename(columns={'j':'i'})
    totald = pd.merge(left=totalm,right=totalf,how='left',on=['i'])
    totald['yd'] = totald['md']+totald['f']
    totals = va_df[['i','go','m','va']].rename(columns={'go':'ys'})
    total = pd.merge(left=totals,right=totald,how='left',on=['i'])
    total['zd']=total['yd']-total['ys']
    total['zdf']=total['zd']/total['ys']
    print 'Excess demand\n'
    print total[['i','zd','zdf']]
