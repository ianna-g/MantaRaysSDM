#pip install requests

#query the weather api data and get the temp on the results
# get long lat from that day, find closest location and query


import requests
def closest_station(lat,long):
    url = f"https://api.weather.gov/points/{lat},{long}"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()
    stations_url = data['properties']['observationStations']

    #List of stations
    stations_response = requests.get(stations_url)
    stations_response.raise_for_status()
    stations_data = stations_response.json()
    #this gets the first station ID and it should be the closest
    station_id = stations_data['features'][0]['properties']['stationIdentifier']
    return station_id

#this is the first lat/long from the 2020 file. later one I'll access the file to do this automatically
#for now I'm just testing to se how API's work
lat = 26.6156
long = -80.0491
date = "2025- 10-06"

try:
    station = closest_station(lat, long)
    print("The nearest station is " + station)

except Exception as e:
    print("Somethign went wrong")
#save the results