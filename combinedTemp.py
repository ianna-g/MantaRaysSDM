#pip install requests
#pip install pandas

#query the weather api data and get the temp on the results
# get long lat from that day, find closest location and query

#import pandas as pd 

#dataBase = pd.read_csv("TrimmedAerial2020.csv")
import requests
def closest_station(lat,long):
    url = f"https://api.weather.gov/points/{lat},{long}"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()
    # gets all the weather observation stations near the long lat that i unputted
    stations_url = data['properties']['observationStations']

    #List of stations
    stations_response = requests.get(stations_url)
    stations_response.raise_for_status()
    stations_data = stations_response.json()
    #this gets the first station ID and it should be the closest
    station_id = stations_data['features'][0]['properties']['stationIdentifier']
    return station_id




#to get the temperature from a specific date

from datetime import datetime, timedelta
def get_temp(station_id, date):
    start = f"{date}T00:00:00Z"
    end = f"{date}T23:59:59Z"
    url = f"https://api.weather.gov/stations/{station_id}/observations?start={start}&end={end}"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    observations = data.get('features', [])
    if not observations:
        print("there was nothing found for this day")
        return None
    temps = []
    for obs in observations:
        temp = obs['properties'].get('temperature')
        if temp and temp['value'] is not None:
            temps.append(temp['value'])

    if not temps:
        print("No temp found")
        return None
    
    averageTemp = sum(temps)/len(temps)
    return averageTemp

#this is the first lat/long from the 2020 file. later one I'll access the file to do this automatically
#for now I'm just testing to se how API's work
lat = 26.6156
long = -80.0491
date = "2025-10-03"

try:
    station = closest_station(lat, long)
    print(f"The nearest station is {station}" )
    temps = get_temp(station, date)
    print(f"the temperature was {temps}°C")

except Exception as e:
    print("Something went wrong")


#save the results