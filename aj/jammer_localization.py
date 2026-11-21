"""Jammer direction finding and localization with the 7+1 CRPA and the INS.

Electronic-support analysis on the stored raw 8-element IF of a campaign item
(<AJ_DATA_ROOT>/campaign/rawIF/<id>) and the exported truth
(<AJ_DATA_ROOT>/campaign/loc/loc_<id>.mat, aj/export_localization.m).
AJ_DATA_ROOT is an environment variable; if unset, the parent folder of the
repository is used (as aj/ajPath.m).

Per block (1 ms of samples every 100 ms while a jammer is on):
  1. sample covariance R of the 8 elements;
  2. source enumeration from the eigenvalues (above 20x the noise floor and
     within 30 dB of the largest, 1 to 6 sources);
  3. directions of the local NED lower hemisphere (ground emitters, elevation
     -45 to 0 deg), rotated to the body frame with the INS attitude (truth +
     attitude error), with the calibrated array model (element positions and
     coupling): signal-subspace scan on a 1 x 0.5 deg grid, then Nelder-Mead
     refinement of the MUSIC null spectrum; mirror candidates kept;
  4. nearest-neighbour association of the bearings into tracks.
Each track is then localized as a stationary ground emitter (height 0 m, a
known-altitude constraint) by Gauss-Newton on azimuth and elevation from the
aircraft positions along the flight path; the estimate is recomputed every
second with all bearings so far. Attitude-error cases show how INS attitude
errors bias the bearings and the location.

Usage:  python -I jammer_localization.py S1_baseline_cruise_s01 [more ids]
Writes  <AJ_DATA_ROOT>/campaign/loc/locres_<id>.json
"""
import json
import os
import sys

import numpy as np
import scipy.io as sio
from scipy.optimize import minimize

FS = 5e6
NMS = 5000                       # samples per ms
BASE = os.environ.get('AJ_DATA_ROOT') or os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(BASE, 'campaign', 'rawIF')
LOC = os.path.join(BASE, 'campaign', 'loc')
A_WGS, F_WGS = 6378137.0, 1 / 298.257223563
E2 = F_WGS * (2 - F_WGS)


def cnb(roll, pitch, head):
    """NED -> body for Euler angles in degrees (as ArrayProcessor/arrayProcess.m)."""
    sp, cp = np.sin(np.radians(roll)), np.cos(np.radians(roll))
    st, ct = np.sin(np.radians(pitch)), np.cos(np.radians(pitch))
    ss, cs = np.sin(np.radians(head)), np.cos(np.radians(head))
    return np.array([[ct * cs, ct * ss, -st],
                     [-cp * ss + sp * st * cs, cp * cs + sp * st * ss, sp * ct],
                     [sp * ss + cp * st * cs, -sp * cs + cp * st * ss, cp * ct]])


def llh2ecef(lat, lon, h):
    lat, lon = np.radians(lat), np.radians(lon)
    n = A_WGS / np.sqrt(1 - E2 * np.sin(lat) ** 2)
    return np.array([(n + h) * np.cos(lat) * np.cos(lon), (n + h) * np.cos(lat) * np.sin(lon),
                     (n * (1 - E2) + h) * np.sin(lat)])


def ecef2llh(x):
    lon = np.arctan2(x[1], x[0])
    p = np.hypot(x[0], x[1])
    lat = np.arctan2(x[2], p * (1 - E2))
    for _ in range(6):
        n = A_WGS / np.sqrt(1 - E2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - n
        lat = np.arctan2(x[2], p * (1 - E2 * n / (n + h)))
    return np.degrees(lat), np.degrees(lon), h


def cen(lat, lon):
    """ECEF -> NED."""
    sl, cl = np.sin(np.radians(lat)), np.cos(np.radians(lat))
    so, co = np.sin(np.radians(lon)), np.cos(np.radians(lon))
    return np.array([[-sl * co, -sl * so, cl], [-so, co, 0.0], [-cl * co, -cl * so, -sl]])


def azel(u):
    """NED unit vector(s) -> azimuth (deg, from north) and elevation (deg, up positive)."""
    return np.degrees(np.arctan2(u[1], u[0])) % 360, np.degrees(np.arcsin(-u[2]))


def unit_ned(az, el):
    az, el = np.radians(az), np.radians(el)
    return np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), -np.sin(el)])


class Item:
    def __init__(self, iid):
        self.id = iid
        tok = iid.split('_')
        self.profile = '_'.join(tok[1:-1])
        d = sio.loadmat(os.path.join(LOC, f'loc_{iid}.mat'), squeeze_me=True, struct_as_record=False)
        self.t = np.atleast_1d(d['t'])
        self.att = np.atleast_2d(d['att'])
        self.pos = np.atleast_2d(d['posEcef'])
        self.uJamB = np.asarray(d['uJamB']).reshape(3, -1, len(self.t))
        self.posB = d['posB']
        self.C = d['C']
        self.lam = float(d['lambda'])
        self.jamPos = np.atleast_2d(d['jamPos'])
        self.jamType = list(np.atleast_1d(d['jamType']))
        self.jamOn = np.atleast_1d(d['jamOn']).astype(float)
        self.jamOff = np.atleast_1d(d['jamOff']).astype(float)
        el = os.path.join(RAW, iid)
        self.files = [os.path.join(el, f'navic8_e{m}_0IF_Fs5MHz_16bit_IQ_Fighter_{self.profile}_{iid}.bin')
                      for m in range(8)]
        self.mm = [np.memmap(f, dtype=np.int16, mode='r') for f in self.files]

    def block(self, kms, nms=1):
        x = np.empty((8, nms * NMS), complex)
        for m in range(8):
            v = np.asarray(self.mm[m][2 * kms * NMS: 2 * (kms + nms) * NMS], float)
            x[m] = v[0::2] + 1j * v[1::2]
        return x

    def att_at(self, tsec):
        """Euler angles (deg) interpolated to tsec; roll and heading unwrapped
        (the 10-ms truth records are too coarse at roll rates of several 100 deg/s)."""
        if not hasattr(self, '_attu'):
            self._attu = np.degrees(np.unwrap(np.radians(self.att), axis=0))
        return np.array([np.interp(tsec, self.t, self._attu[:, i]) for i in range(3)])

    def pos_at(self, tsec):
        return np.array([np.interp(tsec, self.t, self.pos[:, i]) for i in range(3)])

    def rec(self, tsec):
        return int(np.clip(np.searchsorted(self.t, tsec), 0, len(self.t) - 1))

    def steer(self, uB):
        return self.C @ np.exp(-1j * 2 * np.pi / self.lam * (self.posB.T @ uB))


def music_peaks(item, R, Cnb_est, nmax=6):
    w, V = np.linalg.eigh(R)                          # ascending
    sig2 = w[0]
    # sources: eigenvalues above 20x the noise floor and within 30 dB of the
    # strongest (a broadband jammer at very high JNR leaks into secondary
    # eigenvalues through its decorrelation across the aperture)
    D = int(np.clip(np.sum((w > 20 * sig2) & (w > 1e-3 * w[-1])), 1, min(nmax, 7)))
    En = V[:, :8 - D]
    Es = V[:, 8 - D:]

    def steer_n(az, el):
        a = item.steer(Cnb_est @ unit_ned(az, el))
        return a / np.linalg.norm(a, axis=0)

    def spec_s(az, el):
        # signal-subspace projection (smooth, broad peaks): coarse search
        return np.sum(np.abs(Es.conj().T @ steer_n(az, el)) ** 2, axis=0)

    def null_f(x):
        # MUSIC null spectrum (very sharp at high JNR): refinement
        return float(np.sum(np.abs(En.conj().T @ steer_n(np.array([x[0]]), np.array([x[1]]))) ** 2))

    AZ, EL = np.meshgrid(np.arange(0, 360, 1.0), np.arange(-45, 0.01, 0.5))
    P = spec_s(AZ.ravel(), EL.ravel()).reshape(AZ.shape)
    # local maxima (azimuth wraps)
    nb = [np.roll(np.roll(P, i, 0), j, 1) for i in (-1, 0, 1) for j in (-1, 0, 1) if (i, j) != (0, 0)]
    pk = np.all([P >= q for q in nb], axis=0)
    pk[0, :] &= P[0, :] >= P[1, :]
    pk[-1, :] &= P[-1, :] >= P[-2, :]
    idx = np.argsort(P[pk])[::-1]
    cand = list(zip(AZ[pk][idx], EL[pk][idx]))

    def mirror(a, e):
        # the planar array cannot tell a direction from its mirror image in the
        # array (body x-y) plane: both have the same steering vector
        uB = Cnb_est @ unit_ned(a, e)
        uB[2] = -uB[2]
        return azel(Cnb_est.T @ uB)

    def near(a1, e1, a2, e2, tol):
        return min(abs(a1 - a2) % 360, 360 - abs(a1 - a2) % 360) < tol and abs(e1 - e2) < tol

    sets = []                                         # each: list of 1-2 ambiguous NED directions
    for az0, el0 in cand:
        if len(sets) >= D:
            break
        if any(near(az0, el0, a, e, 3) for s in sets for a, e in s):
            continue
        r = minimize(null_f, [az0, el0], method='Nelder-Mead',
                     options={'xatol': 1e-5, 'fatol': 1e-14, 'maxiter': 400, 'initial_simplex':
                              [[az0, el0], [az0 + 0.5, el0], [az0, el0 + 0.5]]})
        a, e = r.x[0] % 360, r.x[1]
        am, em = mirror(a, e)
        alts = [(a, e)] + ([(am, em)] if not near(a, e, am, em, 1.0) else [])
        # directions on the horizon edge of the search region are boundary
        # artefacts (elevation is poorly resolved near the array plane)
        alts = [(x, y) for x, y in alts if -45 < y < -0.3]
        if alts:
            sets.append(alts)
    return sets, D, w / sig2


def ned_at(x):
    lat, lon, _ = ecef2llh(x)
    return cen(lat, lon)


def localize(meas, h_jam=0.0, init=None, sig=(0.05, 0.2)):
    """Gauss-Newton for a ground emitter at height h_jam from (aircraft ECEF, az, el) bearings."""
    P = np.array([m[0] for m in meas])
    za = np.array([m[1] for m in meas])
    ze = np.array([m[2] for m in meas])
    lat0, lon0 = init
    W = np.r_[np.full(len(za), 1 / sig[0] ** 2), np.full(len(ze), 1 / sig[1] ** 2)]
    Ce = [ned_at(p) for p in P]

    def model(lat, lon):
        x = llh2ecef(lat, lon, h_jam)
        a, e = [], []
        for p, C in zip(P, Ce):
            u = C @ (x - p)
            u /= np.linalg.norm(u)
            ai, ei = azel(u)
            a.append(ai)
            e.append(ei)
        return np.array(a), np.array(e)

    lat, lon = lat0, lon0
    mpd_lat = 111320.0
    for _ in range(15):
        a, e = model(lat, lon)
        r = np.r_[(za - a + 180) % 360 - 180, ze - e]
        mpd_lon = 111320.0 * np.cos(np.radians(lat))
        J = np.zeros((len(r), 2))
        for k, (dl, dn) in enumerate(((1 / mpd_lat, 0), (0, 1 / mpd_lon))):  # 1-m steps north, east
            a2, e2 = model(lat + dl, lon + dn)
            J[:, k] = np.r_[(a2 - a + 180) % 360 - 180, e2 - e]
        JW = J * W[:, None]
        Pcov = np.linalg.inv(J.T @ JW)
        dx = Pcov @ (JW.T @ r)
        lat += dx[0] / mpd_lat
        lon += dx[1] / mpd_lon
        if np.hypot(*dx) < 0.01:
            break
    a, e = model(lat, lon)
    res = np.r_[(za - a + 180) % 360 - 180, ze - e]
    return lat, lon, Pcov, res


def cep50(Pcov):
    s = np.sqrt(np.clip(np.linalg.eigvalsh(Pcov), 0, None))
    return 0.59 * (s[0] + s[1]) if s.min() > 0.3 * s.max() else 0.67 * s.max() + 0.23 * s.min()


def run(iid, att_errs):
    it = Item(iid)
    on = it.jamOn[np.isfinite(it.jamOn)]
    t0 = max(15.0, float(np.max(np.r_[on, 0]))) if np.any(it.jamOn > 0) else 0.0
    t1 = float(np.min(np.r_[it.jamOff[np.isfinite(it.jamOff)], 59.0]))
    tb = np.arange(t0 + 0.05, t1 - 0.05, 0.1)
    covs = []
    for tt in tb:
        X = it.block(int(tt * 1000))
        covs.append(X @ X.conj().T / X.shape[1])
    # truth bearings (NED at the aircraft) of every jammer at the block times
    jamE = [llh2ecef(*p) for p in it.jamPos]
    out = {'id': iid, 'jammers': [{'type': str(ty), 'lat': float(p[0]), 'lon': float(p[1])}
                                  for ty, p in zip(it.jamType, it.jamPos)], 'cases': {}}
    for name, err in att_errs.items():
        Cerr = cnb(*err)
        tracks = []                                   # each: list of (t, pos, az, el, truth bearings)
        nsrc = []
        for tt, R in zip(tb, covs):
            tc = tt + 0.5e-3                          # centre of the 1-ms block
            Cnb_est = Cerr @ cnb(*it.att_at(tc))
            sets, D, ev = music_peaks(it, R, Cnb_est)
            nsrc.append(D)
            p = it.pos_at(tc)
            Ce = ned_at(p)
            tru = []
            for xj in jamE:
                u = Ce @ (xj - p)
                tru.append(azel(u / np.linalg.norm(u)))
            used = set()
            for trk in tracks:                        # nearest-neighbour association over all
                if trk[-1][0] < tt - 1.0:             # ambiguous alternatives; tracks unseen
                    continue                          # for 1 s are closed
                la, le = trk[-1][2], trk[-1][3]
                best, bd = None, 3.0
                for i, alts in enumerate(sets):
                    if i in used:
                        continue
                    for a, e in alts:
                        dd = min(abs(a - la) % 360, 360 - abs(a - la) % 360) + 0.5 * abs(e - le)
                        if dd < bd:
                            best, bd = (i, a, e), dd
                if best is not None:
                    used.add(best[0])
                    trk.append((tt, p, best[1], best[2], tru))
            for i, alts in enumerate(sets):           # a new source opens one track per alternative;
                if i not in used:                     # only the true one meets a fixed ground point
                    for a, e in alts:
                        tracks.append([(tt, p, a, e, tru)])
        case = {'att_err_deg': list(err), 'n_blocks': len(tb), 'n_sources_median': float(np.median(nsrc)),
                'tracks': []}
        for trk in tracks:
            if len(trk) < 30:
                continue
            # initial guess: range from the mean elevation (ground emitter), else 30 km
            a0, e0, p0 = trk[0][2], trk[0][3], trk[0][1]
            lat_a, lon_a, h_a = ecef2llh(p0)
            rng0 = h_a / np.tan(np.radians(max(-e0, 2.0)))
            dn, de_ = rng0 * np.cos(np.radians(a0)), rng0 * np.sin(np.radians(a0))
            init = (lat_a + dn / 111320.0, lon_a + de_ / (111320.0 * np.cos(np.radians(lat_a))))
            hist = []
            est = init
            for kend in range(10, len(trk) + 1, 10):  # every second
                meas = [(m[1], m[2], m[3]) for m in trk[:kend]]
                lat, lon, Pc, res = localize(meas, init=est)
                est = (lat, lon)
                x = llh2ecef(lat, lon, 0.0)
                errs = [float(np.linalg.norm(x - xj)) for xj in jamE]
                j = int(np.argmin(errs))
                hist.append({'t': float(trk[kend - 1][0]), 'err_m': errs[j], 'jammer': j, 'cep50_m': float(cep50(Pc))})
            sa = float(np.std(res[:len(res) // 2]))
            se = float(np.std(res[len(res) // 2:]))
            # final pass weighted with the empirical residual spread
            meas = [(m[1], m[2], m[3]) for m in trk]
            lat, lon, Pc, res = localize(meas, init=est, sig=(max(sa, 1e-3), max(se, 1e-3)))
            x = llh2ecef(lat, lon, 0.0)
            errs = [float(np.linalg.norm(x - xj)) for xj in jamE]
            j = int(np.argmin(errs))
            under100 = next((h['t'] for h in hist if h['err_m'] < 100), None)
            case['tracks'].append({'n': len(trk), 't0': float(trk[0][0]), 't1': float(trk[-1][0]),
                                   'az0': float(a0), 'el0': float(e0), 'jammer': j, 'final_err_m': errs[j],
                                   'cep50_m': float(cep50(Pc)), 'res_sd_az_deg': sa, 'res_sd_el_deg': se,
                                   't_first_below_100m': under100, 'history': hist,
                                   'lat': float(lat), 'lon': float(lon), '_meas': trk, '_meas_el': trk})
        # keep the tracks whose bearings meet at one ground point (residual spread
        # in azimuth below 0.2 deg); of kept tracks within 2 km keep the longest
        keep = sorted([c for c in case['tracks'] if c['res_sd_az_deg'] < 0.2], key=lambda c: -c['n'])
        sel = []
        for c in keep:
            if all(np.linalg.norm(llh2ecef(c['lat'], c['lon'], 0) - llh2ecef(q['lat'], q['lon'], 0)) > 2000 for q in sel):
                sel.append(c)
        for c in case['tracks']:
            c['kept'] = any(c is q for q in sel)
        case['doa_rms_az_deg'] = {}
        case['doa_rms_el_deg'] = {}
        for c in sel:
            j = c['jammer']
            ea = np.array([(m[2] - m[4][j][0] + 180) % 360 - 180 for m in c.pop('_meas')])
            ee = np.array([m[3] - m[4][j][1] for m in c['_meas_el']])
            case['doa_rms_az_deg'][j] = float(np.sqrt(np.mean(ea ** 2)))
            case['doa_rms_el_deg'][j] = float(np.sqrt(np.mean(ee ** 2)))
        for c in case['tracks']:
            c.pop('_meas', None);  c.pop('_meas_el', None)
        out['cases'][name] = case
        tr = [f"{c['jammer']}:{c['final_err_m']:.0f}m(n {c['n']}, sd {c['res_sd_az_deg']:.3f})" for c in sel]
        print(f"{iid} {name:8s} D~{case['n_sources_median']:.0f} kept {tr} | DOA rms az "
              f"{ {k: round(v, 4) for k, v in case['doa_rms_az_deg'].items()} } el "
              f"{ {k: round(v, 3) for k, v in case['doa_rms_el_deg'].items()} } | dropped {len(case['tracks']) - len(sel)}", flush=True)
    with open(os.path.join(LOC, f'locres_{iid}.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    cases = {'truth': (0.0, 0.0, 0.0), 'ins': (0.01, -0.01, 0.02), 'yaw0.1': (0.0, 0.0, 0.1), 'yaw0.5': (0.0, 0.0, 0.5)}
    for iid in sys.argv[1:]:
        run(iid, cases)
