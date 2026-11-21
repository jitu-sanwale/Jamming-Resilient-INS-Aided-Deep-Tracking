"""Directions (azimuth, elevation) from the aircraft start over Nashik of the seven-satellite NavIC constellation
used in the simulation: the four real satellites of the deep-tracking study (from broadcast
ephemerides at the start of that evaluation) and three artificial satellites placed at the
vacant nominal slots of the NavIC design (3 GEO + 4 IGSO).

Artificial satellites (nominal geosynchronous geometry, Earth-fixed):
  A1: GEO at 83.0 deg E
  A2, A3: IGSO, inclination 29 deg, equator crossing at 111.75 deg E, arguments of
          latitude 90 and 270 deg (northern and southern lobes of the figure eight).
The GEO model is checked against PRN 6 (GEO 32.5 deg E) and PRN 10 (GEO 129.5 deg E).
"""
import numpy as np

A_WGS, F_WGS = 6378137.0, 1 / 298.257223563
E2 = F_WGS * (2 - F_WGS)
R_GSO = 42164.17e3
SITE = (20.1085, 73.9532, 5000.0)  # aircraft start over Nashik (simulator trajectory): lat, lon (deg), height (m)

REAL = {"PRN 2": (300.5, 66.8), "PRN 6": (246.7, 37.1), "PRN 9": (203.7, 30.8), "PRN 10": (101.2, 24.8)}


def ecef_site(lat, lon, h):
    la, lo = np.radians(lat), np.radians(lon)
    n = A_WGS / np.sqrt(1 - E2 * np.sin(la) ** 2)
    return np.array([(n + h) * np.cos(la) * np.cos(lo), (n + h) * np.cos(la) * np.sin(lo), (n * (1 - E2) + h) * np.sin(la)])


def gso_ecef(incl_deg, lon_cross_deg, u_deg):
    i, u = np.radians(incl_deg), np.radians(u_deg)
    lat = np.arcsin(np.sin(i) * np.sin(u))
    lon = np.radians(lon_cross_deg) + np.arctan2(np.cos(i) * np.sin(u), np.cos(u)) - u
    return R_GSO * np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def azel(sat):
    lat, lon, h = SITE
    r = sat - ecef_site(lat, lon, h)
    la, lo = np.radians(lat), np.radians(lon)
    e = np.array([-np.sin(lo), np.cos(lo), 0.0])
    n = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    u = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    re, rn, ru = r @ e, r @ n, r @ u
    return np.degrees(np.arctan2(re, rn)) % 360, np.degrees(np.arctan2(ru, np.hypot(re, rn)))


ARTIFICIAL_DEF = {"A1 (GEO 83.0E)": (0.0, 83.0, 0.0),
                  "A2 (IGSO 111.75E, N)": (29.0, 111.75, 90.0),
                  "A3 (IGSO 111.75E, S)": (29.0, 111.75, 270.0)}
ARTIFICIAL = {k: tuple(round(float(x), 1) for x in azel(gso_ecef(*v))) for k, v in ARTIFICIAL_DEF.items()}

if __name__ == "__main__":
    for name, lon in (("check PRN 6 (GEO 32.5E)", 32.5), ("check PRN 10 (GEO 129.5E)", 129.5)):
        print(name, [round(float(x), 1) for x in azel(gso_ecef(0.0, lon, 0.0))])
    for k, v in ARTIFICIAL.items():
        print(k, v)
