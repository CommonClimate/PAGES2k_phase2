# Python code to assign country name to all records 
# uses geopy: https://github.com/geopy/geopy
# pandas and scipy (who doesn't?)

import numpy as np
import scipy.io as io
import pandas as pd
import time
from geopy.geocoders import Nominatim, GoogleV3
from geopy.exc import GeocoderTimedOut


geoNom = Nominatim()
#geoogle = GoogleV3()
# load database and extract arrays
P = io.loadmat('../data/PAGES2k_v2.0.0_unpack.mat')
lat = np.squeeze(P['p_lat'])
lon = np.squeeze(P['p_lon'])
names = np.squeeze(P['recordNames'])
nr = len(lon) 

country = ["?" for x in range(nr)]

for r in np.arange(nr):  
    print('Locating rec# ' + str(r+1) + '/' + str(nr) + ', ' + names[r][0])
    loc = geoNom.reverse((lat[r],lon[r]),language='en',timeout=5) # space out requests to avoid timeouts https://operations.osmfoundation.org/policies/nominatim/
    if loc.address: # if geocoding was successful
        add = loc.raw['address']
        if 'country' in add.keys():
            country[r] = add['country']
        elif 'continent' in add.keys():
            country[r] = add['continent'] # mostly for Antartica
        else:
            print("cannot geolocate this record")
    #address = geoogle.reverse((lat[r],lon[r]))

with open('countryNames.csv', 'w') as a_file:
    for result in country:
        result = ''.join(result)
        a_file.write(result + '\n')


# the code below is powerful, but spits out the names in columns
#import xlsxwriter
#
#workbook = xlsxwriter.Workbook('country_names.xlsx')
#worksheet = workbook.add_worksheet()
#row = 0
#
#for col, data in enumerate(np.array(country)):
#    worksheet.write_column(row, col, data)
#
#workbook.close()

