pip install geopy
from geopy.geocoders import Nominatim
geolocator = Nominatim()
location = geolocator.reverse("48.8588443, 2.2943506")

print(location.address)
print (location.raw)
