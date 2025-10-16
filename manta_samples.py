#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 09:34:31 2025

@author: EcoVision Analytics
"""

"""
Performs per-row extraction based on : Trimmed_Aerial_Data_backup_file_08302025.csv
  • Distance to nearest inlet (km)
  • NOAA daily air temperature (sample day)
  • NOAA precipitation (previous day)
  • Lunar phase + illumination
  • NOAA CO-OPS northbound current velocity (m/s; +northbound, -southbound), daily mean
      - Observations first; if none, falls back to CO-OPS current *predictions*
  • Tide information at sampling time from CO-OPS tide predictions:
      - predicted tide height (meters)
      - tide stage: Flood (Incoming), Ebb (Outgoing), or Slack
      - intensity: weak / moderate / strong (based on slope magnitude)

INPUT (required columns): Latitude, Longitude, Time  (Time e.g. 2020-01-08T15:05:11Z)

Usage:
  1) pip install pandas numpy requests
  2) python mmf_samples.py
"""

# --------------------------- User Config ------------------------------------
# Assumes PC operations. Adjust accordingly for Mac. Optionally set your working directory (uncomment and edit):
import os
os.chdir(r"Set\your\CWD")

IN_CSV  = r"Trimmed_Aerial_Data_backup_file_08302025.csv" # Renamed
OUT_CSV = r"Trimmed_Aerial_Data_draft.csv"

# NOAA CDO token. Obtain token here: https://www.ncdc.noaa.gov/cdo-web/token enter email, check email, input below
NOAA_TOKEN = "YoUrToKeNhErE"
# ---------------------------------------------------------------------------

BIN_DEG = 0.1

# Max acceptable CO-OPS station distance (km)
MAX_COOPS_DIST_KM = 80.0

# Parallel workers for (station,date) fetches.
# Keep modest to avoid rate limits; set to 1 to disable threading.
WORKERS_GHCND = 4
WORKERS_COOPS = 4

VERBOSE = True

# --------------------------- Imports ----------------------------------------
import math
import time
import sys
from time import monotonic
from datetime import datetime, timedelta, timezone
from typing import Dict, Optional, Tuple, List, Iterable

import pandas as pd
import numpy as np
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from concurrent.futures import ThreadPoolExecutor, as_completed

# --------------------------- Logging ----------------------------------------
def log(msg: str):
    if VERBOSE:
        print(msg)
        sys.stdout.flush()

# --------------------------- NOAA endpoints ---------------------------------
CDO_BASE    = "https://www.ncdc.noaa.gov/cdo-web/api/v2"
COOPS_MDAPI = "https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations.json"
COOPS_DATA  = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"

# Timings / retries
CDO_SLEEP = 0.05
COOPS_SLEEP = 0.05
CDO_TIMEOUT = (5, 45)
COOPS_TIMEOUT = (5, 45)

# --------------------------- HTTP session -----------------------------------
def _make_session() -> requests.Session:
    s = requests.Session()
    retry = Retry(
        total=4, connect=4, read=4,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET"]),
        raise_on_status=False,
        respect_retry_after_header=True,
    )
    ad = HTTPAdapter(max_retries=retry, pool_connections=40, pool_maxsize=40)
    s.mount("https://", ad)
    s.mount("http://", ad)
    return s

_HTTP = _make_session()

def _cdo_get(path: str, params: Dict) -> Dict:
    url = f"{CDO_BASE}/{path}"
    for attempt in range(5):
        try:
            if attempt: time.sleep(CDO_SLEEP * (2 ** (attempt - 1)))
            r = _HTTP.get(url, headers={"token": NOAA_TOKEN}, params=params, timeout=CDO_TIMEOUT)
            r.raise_for_status()
            return r.json()
        except requests.exceptions.RequestException as e:
            err = e
            continue
    raise err

def _coops_get(url: str, params: Dict) -> Dict:
    for attempt in range(5):
        try:
            if attempt: time.sleep(COOPS_SLEEP * (2 ** (attempt - 1)))
            r = _HTTP.get(url, params=params, timeout=COOPS_TIMEOUT)
            r.raise_for_status()
            return r.json()
        except requests.exceptions.RequestException as e:
            err = e
            continue
    raise err

# --------------------------- Geometry ---------------------------------------
def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0088
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat/2.0)**2 + np.cos(np.radians(lat1))*np.cos(np.radians(lat2))*np.sin(dlon/2.0)**2
    return 2 * R * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

# Vectorized nearest inlet
INLETS = np.array([
    ("Government Cut (Miami)", 25.7660, -80.1325),
    ("Port Everglades Inlet", 26.0945, -80.0907),
    ("Hillsboro Inlet", 26.2533, -80.0763),
    ("Boca Raton Inlet", 26.3583, -80.0803),
    ("Lake Worth Inlet", 26.7708, -80.0300),
    ("Jupiter Inlet", 26.9440, -80.0634),
    ("St. Lucie Inlet", 27.1667, -80.1400),
    ("Fort Pierce Inlet", 27.4750, -80.2700),
    ("Sebastian Inlet", 27.8603, -80.4473),
], dtype=object)

INLET_NAMES = INLETS[:,0]
INLET_LATS  = INLETS[:,1].astype(float)
INLET_LONS  = INLETS[:,2].astype(float)

def nearest_inlet_vec(lats: np.ndarray, lons: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    # compute distances to all inlets, pick min
    dists = np.stack([haversine_km(lats, lons, INLET_LATS[i], INLET_LONS[i]) for i in range(len(INLETS))], axis=1)
    idx = np.argmin(dists, axis=1)
    return INLET_NAMES[idx], dists[np.arange(len(lats)), idx]

# --------------------------- Moon phase (per date) ---------------------------
def moon_phase(date: datetime) -> Tuple[str, float]:
    known_new_moon = datetime(2000, 1, 6, 18, 14, tzinfo=timezone.utc)
    synodic_month = 29.53058867
    days = (date.replace(tzinfo=timezone.utc) - known_new_moon).total_seconds() / 86400.0
    phase = (days / synodic_month) % 1.0
    illum = 0.5 * (1 - math.cos(2 * math.pi * phase))
    if phase < 0.03 or phase > 0.97: name = "New"
    elif phase < 0.22: name = "Waxing Crescent"
    elif phase < 0.28: name = "First Quarter"
    elif phase < 0.47: name = "Waxing Gibbous"
    elif phase < 0.53: name = "Full"
    elif phase < 0.72: name = "Waning Gibbous"
    elif phase < 0.78: name = "Last Quarter"
    else: name = "Waning Crescent"
    return name, float(illum)

# --------------------------- Station lookup & caches -------------------------
def bin_key(lat, lon):
    """Round coordinates to BIN_DEG grid for station caching.
    Returns None if lat/lon are missing or invalid."""
    try:
        if pd.isna(lat) or pd.isna(lon):
            return None
        return (
            round(float(lat) / BIN_DEG) * BIN_DEG,
            round(float(lon) / BIN_DEG) * BIN_DEG
        )
    except Exception:
        return None

_GHCND_STATION_CACHE: Dict[Tuple[float,float], Optional[Dict]] = {}
_COOPS_CURR_STATIONS: Optional[List[Dict]] = None
_COOPS_TIDE_STATIONS: Optional[List[Dict]] = None
_TIDE_DAY_SERIES: Dict[Tuple[str, datetime], Tuple[np.ndarray, np.ndarray]] = {}  # (station, date) -> (tsecs[], vals[])
_GHCND_DAY_CACHE: Dict[Tuple[str, datetime], Tuple[Optional[float], Optional[float]]] = {}
_COOPS_CURR_DAY_CACHE: Dict[Tuple[str, datetime], Optional[float]] = {}

def get_coops_stations_of_type(stype: str) -> List[Dict]:
    js = _coops_get(COOPS_MDAPI, {"type": stype})
    return js.get("stations", []) or js.get("results", []) or []

def nearest_coops_station(lat: float, lon: float, stype: str) -> Optional[Dict]:
    global _COOPS_CURR_STATIONS, _COOPS_TIDE_STATIONS
    table = None
    if stype == "currents":
        if _COOPS_CURR_STATIONS is None:
            _COOPS_CURR_STATIONS = get_coops_stations_of_type("currents")
        table = _COOPS_CURR_STATIONS
    else:
        if _COOPS_TIDE_STATIONS is None:
            _COOPS_TIDE_STATIONS = get_coops_stations_of_type("waterlevels")
        table = _COOPS_TIDE_STATIONS
    best, best_d = None, float("inf")
    for s in table:
        try:
            slat = float(s.get("lat", s.get("latitude")))
            slon = float(s.get("lng", s.get("longitude")))
        except Exception:
            continue
        d = haversine_km(lat, lon, slat, slon)
        if d < best_d:
            best, best_d = s, d
    if best and best_d <= MAX_COOPS_DIST_KM:
        best["_distance_km"] = best_d
        return best
    return None

def find_nearest_ghcnd_station(lat: float, lon: float, bbox_size: float = 0.5) -> Optional[Dict]:
    key = bin_key(lat, lon)
    if key in _GHCND_STATION_CACHE:
        return _GHCND_STATION_CACHE[key]
    best = None
    for mult in [1, 2, 4, 8]:
        delta = bbox_size * mult
        extent = f"{lat - delta},{lon - delta},{lat + delta},{lon + delta}"
        js = _cdo_get("stations", {"datasetid": "GHCND", "extent": extent, "limit": 1000})
        stations = js.get("results", [])
        if not stations: continue
        best = min(stations, key=lambda s: haversine_km(lat, lon, s["latitude"], s["longitude"]))
        break
    _GHCND_STATION_CACHE[key] = best
    return best

# --------------------------- GHCND per (station,date) ------------------------
def ghcnd_daily(station_id: str, date: datetime) -> Tuple[Optional[float], Optional[float]]:
    day_key = (station_id, date)
    if day_key in _GHCND_DAY_CACHE:
        return _GHCND_DAY_CACHE[day_key]
    day = date.strftime("%Y-%m-%d")
    prev = (date - timedelta(days=1)).strftime("%Y-%m-%d")
    tavg_c, prcp_mm = None, None
    try:
        vals = _cdo_get("data", {
            "datasetid": "GHCND",
            "datatypeid": ["TAVG","TMAX","TMIN"],
            "stations": station_id, "startdate": day, "enddate": day,
            "units": "metric", "limit": 1000
        }).get("results", [])
        d = {v.get("datatype"): v.get("value") for v in vals}
        if d.get("TAVG") is not None:
            tavg_c = float(d["TAVG"])
        elif d.get("TMAX") is not None and d.get("TMIN") is not None:
            tavg_c = (float(d["TMAX"]) + float(d["TMIN"])) / 2.0
    except Exception:
        pass
    try:
        vals = _cdo_get("data", {
            "datasetid": "GHCND",
            "datatypeid": "PRCP",
            "stations": station_id, "startdate": prev, "enddate": prev,
            "units": "metric", "limit": 1000
        }).get("results", [])
        if vals:
            prcp_mm = float(vals[0].get("value"))
    except Exception:
        pass
    _GHCND_DAY_CACHE[day_key] = (tavg_c, prcp_mm)
    return tavg_c, prcp_mm

# --------------------------- Currents per (station,date) ---------------------
def coops_currents_obs_daily_mean_north(station_id: str, date: datetime) -> Optional[float]:
    p = {
        "product": "currents",
        "application": "mmf-enrichment",
        "begin_date": date.strftime("%Y%m%d 00:00"),
        "end_date":   date.strftime("%Y%m%d 23:59"),
        "time_zone": "gmt", "interval": "1",
        "units": "metric", "format": "json",
        "station": station_id,
    }
    try:
        js = _coops_get(COOPS_DATA, p)
        data = js.get("data", [])
        if not data: return None
        comps = []
        for r in data:
            try:
                v = float(r["Speed"])
                deg = float(r["Direction"])
                comps.append(v * math.cos(math.radians(deg)))
            except Exception:
                pass
        return float(np.nanmean(comps)) if comps else None
    except Exception:
        return None

def coops_currents_pred_daily_mean_north(station_id: str, date: datetime) -> Optional[float]:
    p = {
        "product": "currents_predictions",
        "application": "mmf-enrichment",
        "begin_date": date.strftime("%Y%m%d 00:00"),
        "end_date":   date.strftime("%Y%m%d 23:59"),
        "time_zone": "gmt", "interval": "6",
        "units": "metric", "format": "json",
        "station": station_id,
    }
    try:
        js = _coops_get(COOPS_DATA, p)
        data = js.get("data", [])
        if not data: return None
        comps = []
        for r in data:
            try:
                v = float(r["Speed"])
                deg = float(r["Direction"])
                comps.append(v * math.cos(math.radians(deg)))
            except Exception:
                pass
        return float(np.nanmean(comps)) if comps else None
    except Exception:
        return None

def daily_mean_north_component(station_id: str, date: datetime) -> Optional[float]:
    key = (station_id, date)
    if key in _COOPS_CURR_DAY_CACHE:
        return _COOPS_CURR_DAY_CACHE[key]
    north = coops_currents_obs_daily_mean_north(station_id, date)
    if north is None:
        north = coops_currents_pred_daily_mean_north(station_id, date)
    _COOPS_CURR_DAY_CACHE[key] = north
    return north

# --------------------------- Tide series per (station,date) ------------------
def fetch_tide_day_series(station_id: str, date: datetime) -> Tuple[np.ndarray, np.ndarray]:
    key = (station_id, date)
    if key in _TIDE_DAY_SERIES:
        return _TIDE_DAY_SERIES[key]
    begin = date.strftime("%Y%m%d 00:00")
    end   = date.strftime("%Y%m%d 23:59")
    p = {
        "product": "predictions",
        "application": "mmf-enrichment",
        "begin_date": begin,
        "end_date": end,
        "time_zone": "gmt",
        "interval": "6",
        "units": "metric",
        "format": "json",
        "station": station_id,
    }
    try:
        js = _coops_get(COOPS_DATA, p)
        preds = js.get("predictions", []) or js.get("data", [])
        if not preds: 
            arr_t = np.array([])
            arr_v = np.array([])
        else:
            arr_t = np.array([datetime.strptime(p["t"], "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc).timestamp() for p in preds], dtype=float)
            arr_v = np.array([float(p["v"]) for p in preds], dtype=float)
    except Exception:
        arr_t = np.array([])
        arr_v = np.array([])
    _TIDE_DAY_SERIES[key] = (arr_t, arr_v)
    return arr_t, arr_v

def tide_at_time_from_series(arr_t: np.ndarray, arr_v: np.ndarray, dt: datetime) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    if arr_t.size < 3: 
        return None, None, None
    ts = dt.timestamp()
    idx = np.searchsorted(arr_t, ts)
    if idx == 0:
        h = arr_v[0]
    elif idx >= arr_t.size:
        h = arr_v[-1]
    else:
        t0, t1 = arr_t[idx-1], arr_t[idx]
        v0, v1 = arr_v[idx-1], arr_v[idx]
        w = (ts - t0) / max((t1 - t0), 1.0)
        h = v0 + w * (v1 - v0)
    # slope using +/- 6 minutes (1 sample step)
    i0 = max(0, idx-1)
    i1 = min(arr_t.size-1, idx)
    # widen a bit for robustness
    i0 = max(0, i0-1)
    i1 = min(arr_t.size-1, i1+1)
    dt_hr = max((arr_t[i1] - arr_t[i0]) / 3600.0, 1e-6)
    slope = (arr_v[i1] - arr_v[i0]) / dt_hr  # m/hour
    eps = 0.005
    if slope > eps: stage = "Flood (Incoming)"
    elif slope < -eps: stage = "Ebb (Outgoing)"
    else: stage = "Slack"
    a = abs(slope)
    intensity = "weak" if a < 0.05 else ("moderate" if a < 0.15 else "strong")
    return float(h), stage, intensity

# --------------------------- Batch helpers -----------------------------------
def unique_dates(datetimes: pd.Series) -> List[datetime]:
    dates = datetimes.dt.floor("D").dropna().unique()
    return [datetime(d.year, d.month, d.day, tzinfo=timezone.utc) for d in pd.to_datetime(dates)]

def group_by_bin(df: pd.DataFrame) -> pd.Series:
    return df.apply(lambda r: bin_key(r["Latitude"], r["Longitude"]), axis=1)

# ------------------------------- Main ----------------------------------------
def main():
    t0 = monotonic()
    df = pd.read_csv(IN_CSV)
    df.columns = [c.strip() for c in df.columns]
    for required in ("Latitude", "Longitude", "Time"):
        if required not in df.columns:
            raise ValueError(f"Missing required column: {required}")

    # Parse datetimes
    df["_dt"] = pd.to_datetime(df["Time"], errors="coerce")
    if df["_dt"].isna().any():
        log("WARNING: Some 'Time' values could not be parsed; those rows may produce NaN outputs.")
    if df["_dt"].dt.tz is None:
        df["_date_utc"] = df["_dt"].dt.tz_localize("UTC").dt.floor("D")
    else:
        df["_date_utc"] = df["_dt"].dt.tz_convert("UTC").dt.floor("D") 

    # Vectorized nearest inlet
    names, dkm = nearest_inlet_vec(df["Latitude"].to_numpy(float), df["Longitude"].to_numpy(float))

    # Prepare output columns with proper dtypes
    df["nearest_inlet_name"] = pd.Series(names, dtype="string")
    df["nearest_inlet_km"]   = pd.Series(dkm, dtype="float")
    df["lunar_phase"]        = pd.Series(dtype="string")
    df["lunar_illum"]        = pd.Series(dtype="float")
    df["tavg_c"]             = pd.Series(dtype="float")
    df["prcp_mm_prev_day"]   = pd.Series(dtype="float")
    df["curr_north_mps"]     = pd.Series(dtype="float")
    df["tide_height_m"]      = pd.Series(dtype="float")
    df["tide_stage"]         = pd.Series(dtype="string")
    df["tide_intensity"]     = pd.Series(dtype="string")

    # --------------------- Precompute lunar phase per unique date -------------
    dates = unique_dates(df["_dt"])
    moon_map: Dict[datetime, Tuple[str,float]] = {d: moon_phase(d) for d in dates}
    df["lunar_phase"] = df["_date_utc"].map(lambda d: None if pd.isna(d) else moon_map[datetime(d.year,d.month,d.day,tzinfo=timezone.utc)][0])
    df["lunar_illum"] = df["_date_utc"].map(lambda d: np.nan if pd.isna(d) else moon_map[datetime(d.year,d.month,d.day,tzinfo=timezone.utc)][1])

    # --------------------- Resolve stations per lat/lon bin -------------------
    log("Resolving nearest stations per lat/lon bin…")
    df["_bin"] = group_by_bin(df)
    bin_to_ghcnd = {}
    bin_to_curr  = {}
    bin_to_tide  = {}
    for b in df["_bin"].dropna().unique():
        lat_b, lon_b = b
        st_g = find_nearest_ghcnd_station(lat_b, lon_b)
        bin_to_ghcnd[b] = None if not st_g else st_g.get("id")
        st_c = nearest_coops_station(lat_b, lon_b, "currents")
        bin_to_curr[b]  = None if not st_c else (st_c.get("id") or st_c.get("station"))
        st_t = nearest_coops_station(lat_b, lon_b, "waterlevels")
        bin_to_tide[b]  = None if not st_t else (st_t.get("id") or st_t.get("station"))

    df["_ghcnd_id"] = df["_bin"].map(bin_to_ghcnd)
    df["_curr_id"]  = df["_bin"].map(bin_to_curr)
    df["_tide_id"]  = df["_bin"].map(bin_to_tide)

    # --------------------- Batch fetch: GHCND per (station, date) -------------
    ghcnd_tasks = {(sid, d) for sid in df["_ghcnd_id"].dropna().unique()
                         for d in df["_date_utc"].dropna().unique()}
    log(f"GHCND tasks: {len(ghcnd_tasks)}")
    def _do_ghcnd(k):
        sid, d = k
        return k, ghcnd_daily(sid, datetime(d.year,d.month,d.day,tzinfo=timezone.utc))
    if WORKERS_GHCND > 1:
        with ThreadPoolExecutor(max_workers=WORKERS_GHCND) as ex:
            for fut in as_completed([ex.submit(_do_ghcnd, k) for k in ghcnd_tasks]):
                k, val = fut.result()
                _GHCND_DAY_CACHE[k] = val
    else:
        for k in ghcnd_tasks:
            _GHCND_DAY_CACHE[k] = ghcnd_daily(k[0], datetime(k[1].year,k[1].month,k[1].day,tzinfo=timezone.utc))

    # Map back to rows
    def map_ghcnd(row):
        sid = row["_ghcnd_id"]; d = row["_date_utc"]
        if pd.isna(d) or sid is None: return (np.nan, np.nan)
        t, p = _GHCND_DAY_CACHE.get((sid, datetime(d.year,d.month,d.day,tzinfo=timezone.utc)), (None,None))
        return (np.nan if t is None else float(t), np.nan if p is None else float(p))
    tmp = df.apply(map_ghcnd, axis=1, result_type="expand")
    df["tavg_c"] = tmp[0]
    df["prcp_mm_prev_day"] = tmp[1]

    # --------------------- Batch fetch: Currents per (station, date) ----------
    curr_tasks = {(sid, d) for sid in df["_curr_id"].dropna().unique()
                        for d in df["_date_utc"].dropna().unique()}
    log(f"Currents tasks: {len(curr_tasks)}")
    def _do_curr(k):
        sid, d = k
        return k, daily_mean_north_component(sid, datetime(d.year,d.month,d.day,tzinfo=timezone.utc))
    if WORKERS_COOPS > 1:
        with ThreadPoolExecutor(max_workers=WORKERS_COOPS) as ex:
            for fut in as_completed([ex.submit(_do_curr, k) for k in curr_tasks]):
                k, val = fut.result()
                _COOPS_CURR_DAY_CACHE[k] = val
    else:
        for k in curr_tasks:
            _COOPS_CURR_DAY_CACHE[k] = daily_mean_north_component(k[0], datetime(k[1].year,k[1].month,k[1].day,tzinfo=timezone.utc))

    def map_curr(row):
        sid = row["_curr_id"]; d = row["_date_utc"]
        if pd.isna(d) or sid is None: return np.nan
        val = _COOPS_CURR_DAY_CACHE.get((sid, datetime(d.year,d.month,d.day,tzinfo=timezone.utc)))
        return np.nan if val is None else float(val)
    df["curr_north_mps"] = df.apply(map_curr, axis=1)

    # --------------------- Batch fetch: Tide day series per (station, date) ---
    tide_tasks = {(sid, d) for sid in df["_tide_id"].dropna().unique()
                        for d in df["_date_utc"].dropna().unique()}
    log(f"Tide day-series tasks: {len(tide_tasks)}")
    for sid, d in tide_tasks:
        fetch_tide_day_series(sid, datetime(d.year,d.month,d.day,tzinfo=timezone.utc))

    # Interpolate tide per row (fast, local)
    def map_tide(row):
        sid = row["_tide_id"]; dtv = row["_dt"]
        if sid is None or pd.isna(dtv): return (np.nan, None, None)
        d = datetime(dtv.year, dtv.month, dtv.day, tzinfo=timezone.utc)
        arr_t, arr_v = _TIDE_DAY_SERIES.get((sid, d), (np.array([]), np.array([])))
        if arr_t.size == 0: return (np.nan, None, None)
        h, stage, intensity = tide_at_time_from_series(arr_t, arr_v, dtv.tz_localize("UTC") if dtv.tzinfo is None else dtv.tz_convert("UTC"))
        return (np.nan if h is None else float(h), stage, intensity)
    tmp = df.apply(map_tide, axis=1, result_type="expand")
    df["tide_height_m"] = tmp[0]
    df["tide_stage"] = tmp[1].astype("string")
    df["tide_intensity"] = tmp[2].astype("string")

    # --------------------- Cleanup & write ------------------------------------
    df.drop(columns=["_dt","_date_utc","_bin","_ghcnd_id","_curr_id","_tide_id"], inplace=True, errors="ignore")
    df.to_csv(OUT_CSV, index=False)
    log(f" Done in {monotonic()-t0:.1f}s → {OUT_CSV}")

if __name__ == "__main__":
    main()