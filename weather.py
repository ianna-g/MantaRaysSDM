import requests
import math

ERDDAP_BASE = "https://coastwatch.pfeg.noaa.gov/erddap"


def haversine(lat1, lon1, lat2, lon2):
    R = 6371
    d_lat = math.radians(lat2 - lat1)
    d_lon = math.radians(lon2 - lon1)
    a = (math.sin(d_lat/2)**2 +
         math.cos(math.radians(lat1)) *
         math.cos(math.radians(lat2)) *
         math.sin(d_lon/2)**2)
    return 2 * R * math.asin(math.sqrt(a))


def get_all_stations():
    url = f"{ERDDAP_BASE}/tabledap/cwwcNDBCMet.json?station,latitude,longitude"
    r = requests.get(url)
    r.raise_for_status()
    return r.json()["table"]["rows"]


def station_has_temp(station_id):
    """Check if a station has either air or water temperature."""
    test_url = (
        f"{ERDDAP_BASE}/tabledap/cwwcNDBCMet.json?"
        f"air_temperature,sea_surface_temperature&station=\"{station_id}\"&time>=2000-01-01T00:00:00Z&time<=2000-01-02T00:00:00Z"
    )
    r = requests.get(test_url)
    # 400 = station does NOT have variable
    return r.status_code == 200


def closest_valid_station(lat, lon):
    print("Fetching all NDBC stations…")
    stations = get_all_stations()

    print("Filtering to stations that actually have temperature…")
    valid = []

    for stn, s_lat, s_lon in stations:
        if station_has_temp(stn):
            d = haversine(lat, lon, s_lat, s_lon)
            valid.append((d, stn, s_lat, s_lon))

    if not valid:
        print("No temperature-capable stations found.")
        return None

    valid.sort(key=lambda x: x[0])
    d, station_id, s_lat, s_lon = valid[0]

    print(f"Nearest valid temp station: {station_id} ({d:.1f} km)")
    return station_id


def get_air_and_water_temp(station_id, date):
    start = f"{date}T00:00:00Z"
    end   = f"{date}T23:59:59Z"

    url = (
        f"{ERDDAP_BASE}/tabledap/cwwcNDBCMet.json"
        f"?time,air_temperature,sea_surface_temperature"
        f"&station=\"{station_id}\""
        f"&time>={start}"
        f"&time<={end}"
    )

    print(f"\nFetching temperature for station {station_id}…")

    r = requests.get(url)
    r.raise_for_status()
    rows = r.json()["table"]["rows"]

    air = [row[1] for row in rows if row[1] is not None]
    water = [row[2] for row in rows if row[2] is not None]

    return {
        "air": sum(air) / len(air) if air else None,
        "water": sum(water) / len(water) if water else None,
    }


# ---------------------------
# Run
# ---------------------------

lat = 26.6156
lon = -80.0491
date = "2025-11-13"

print("\n=== ERDDAP WEATHER ===")
station = closest_valid_station(lat, lon)

if station:
    temps = get_air_and_water_temp(station, date)
    print("\n=== RESULTS ===")
    print("Station:", station)
    print("Air Temp (°C):", temps["air"])
    print("Water Temp (°C):", temps["water"])
