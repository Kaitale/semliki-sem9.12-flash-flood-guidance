# -*- coding: utf-8 -*-
"""ffg_core.py -- 21 Nov 2026: flash-flood guidance from the HEC-HMS states of Semuliki_Ops_Cal.

The rain (mm in d hours, uniform over the point's subbasins) that brings the point to a target flow, from the state
of each subbasin at the issue: the same Deficit and Constant loss (canopy, surface store, soil deficit, constant
rate, impervious share) and the same Clark unit hydrograph as HEC-HMS (default time-area curve; checked against
the FLOW-UNIT GRAPH HEC-HMS writes: peak 0.813 against 0.817 m3/s per mm for S34), and the percolation routed
through the model's two linear groundwater reservoirs, whose quick layer adds to the peak in a long storm on a
saturated catchment (checked against HEC-HMS with three synthetic storms on the Wassa). The flow at the issue is
added as it stands (conservative: the recession is not subtracted)."""
import re, numpy as np


def parse_basin(text, names):
    out = {}
    for s in names:
        m = re.search(r"(?ms)^Subbasin: %s\n(.*?)^End:" % re.escape(s), text)
        b = m.group(1)
        g = lambda k: float(re.search(r"(?m)^\s+%s: ([0-9.eE+-]+)" % k, b).group(1))
        out[s] = dict(area=g("Area"), imperv=g("Percent Impervious Area") / 100, maxdef=g("Maximum Deficit"), perc=g("Percolation Rate"),
                      tc=g("Time of Concentration"), r=g("Storage Coefficient"), canopy=g("Canopy Storage Capacity"),
                      surface=g("Surface Storage Capacity"), gw1f=g("GW-1 Baseflow Fraction"), gw1k=g("GW-1 Routing Coefficient"),
                      gw2f=g("GW-2 Baseflow Fraction"), gw2k=g("GW-2 Routing Coefficient"))
    return out


def clark_uh(area_km2, tc, r, dt=1.0):
    n = int(tc + 8 * r) + 5
    def A(t):
        x = min(max(t / tc, 0), 1)
        return 1.414 * x ** 1.5 if x <= 0.5 else 1 - 1.414 * (1 - x) ** 1.5
    inflow = np.array([A((i + 1) * dt) - A(i * dt) for i in range(n)]) * area_km2 * 1e3 / (dt * 3600)
    c = dt / (r + 0.5 * dt); o = np.zeros(n); prev = 0.0
    for i in range(n):
        prev = c * inflow[i] + (1 - c) * prev; o[i] = prev
    return np.r_[o[0] / 2, (o[1:] + o[:-1]) / 2]


def excess(rain, p, st, perc_out=None, tail=48):
    """hourly excess (mm) of an hourly rain sequence, from the state st of the subbasin, in the order HEC-HMS uses
    (checked against its PRECIP-EXCESS, PERC-SOIL and STORAGE-SURFACE): canopy interception; the impervious share runs
    off; the soil takes what its deficit holds, then the constant rate; what is left fills the surface store and
    spills as excess; the surface store drains into the soil at the constant rate after the rain. perc_out (a list)
    receives the hourly percolation, which HEC-HMS sends to the two linear-reservoir groundwater layers."""
    can = max(p["canopy"] - st.get("canopy", 0.0), 0.0); sto = min(st.get("surface", 0.0), p["surface"])
    de = st.get("deficit", p["maxdef"]); ex = []
    for fall in list(rain) + [0.0] * tail:
        t = min(can, fall); can -= t; fall -= t
        imp = p["imperv"] * fall; water = fall - imp + sto; sto = 0.0
        t = min(de, water); de -= t; water -= t
        pc = min(water, p["perc"]) if de <= 0 else 0.0; water -= pc
        if perc_out is not None:
            perc_out.append(pc)
        t = min(water, p["surface"]); sto = t; water -= t
        ex.append(imp + water)
    return np.array(ex)


def linres(inflow_mm, area_km2, k_h, n_out):
    """single linear reservoir (HEC-HMS baseflow layer), m3/s from hourly inflow in mm, empty at the start"""
    c = 1.0 / (k_h + 0.5); o = np.zeros(n_out); prev = 0.0
    q_in = np.pad(np.asarray(inflow_mm, float), (0, max(0, n_out - len(inflow_mm)))) * area_km2 * 1e3 / 3600.0
    for i in range(n_out):
        prev = c * q_in[i] + (1 - c) * prev; o[i] = prev
    return o


def response(rain_by_sub, P, S):
    """flow added at the point (m3/s): direct runoff (Clark) plus the quick part of the groundwater response"""
    tot = np.zeros(1)
    for s, p in P.items():
        pc = []; ex = excess(rain_by_sub[s], p, S[s], pc)
        q = np.convolve(ex, clark_uh(p["area"], p["tc"], p["r"]))
        q = q + linres(np.array(pc) * p.get("gw1f", 0.0), p["area"], p.get("gw1k", 1e9), len(q)) \
              + linres(np.array(pc) * p.get("gw2f", 0.0), p["area"], p.get("gw2k", 1e9), len(q))
        n = max(len(tot), len(q)); tot = np.pad(tot, (0, n - len(tot))) + np.pad(q, (0, n - len(q)))
    return tot


def peak(depth, d_h, P, S, base=0.0):
    return base + float(response({s: np.full(int(d_h), depth / d_h) for s in P}, P, S).max())


def guidance(target, d_h, P, S, base, hi=400.0):
    if base >= target:
        return 0.0
    if peak(hi, d_h, P, S, base) < target:
        return None
    lo = 0.0
    for _ in range(40):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if peak(m, d_h, P, S, base) < target else (lo, m)
    return round(hi, 1)
