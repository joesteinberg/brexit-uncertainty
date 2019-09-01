##########################################################################################
# This script contains the regional aggregations used by other scripts in this folder.
# It has no main body, you must import it in the other scripts.
#
# regions is a dictionary, where the key is a region name and the value is a list of countries
# in the region.
#
# regions_dict is another dictionary, where the key is the region name and the value is the
# region number.
#
# The function "which_region" tells you which region a given country belongs to.
#
# The function "which_region_num" tells you the region number associated with a given region
##########################################################################################

regions=[]
regions_dict=[]
path=''

nr=3
path = 'output/'

# UK
ruk=['GBR']

# European Union
reu=['AUT',
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

regions = {'UK':ruk,
           'EU':reu}

regions_dict = {'UK':1,'EU':2,'Constructed RoW':3};

def which_region(country):
    for key,val in regions.iteritems():
        if country in val:
            return key
    return 'ROW'

def which_region_num(region):
    if region in regions:
        return regions_dict[region]
    else:
        return nr

def which_region_name(num):
    for key,val in regions_dict.iteritems():
        if val == num:
            return key
    return ''
