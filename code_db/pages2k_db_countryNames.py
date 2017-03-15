# Python code to assign country name to all records 
# uses geopy: https://github.com/geopy/geopy
# pandas and scipy (who doesn't?)

import numpy as np
import scipy.io as io
import pandas as pd
from geopy.geocoders import Nominatim, GoogleV3
geoNom = Nominatim()
#geoogle = GoogleV3()
# load database and extract arrays
P = io.loadmat('../data/PAGES2k_v2.0.0_unpack.mat')
lat = np.around(np.squeeze(P['p_lat']),4)
lon = np.around(np.squeeze(P['p_lon']),4)
names = np.squeeze(P['recordNames'])
nr = len(lon) 
country = [None] * nr
for r in np.arange(nr):    
    loc = geoNom.reverse((lat[r],lon[r]),language='en') 
    if loc.raw['address']['country']: # if geocoding was successful
        country[r] = loc.raw['address']['country']
    #address = geoogle.reverse((lat[r],lon[r]))

print(country)
