"""Calcul des références chromatographiques sur le point GAM-6.

Le script ne fixe aucune tolérance métier. Il calcule :
- temps de rétention de référence à l'apex du quantifiant ;
- bornes d'intégration quantifiant et qualifiant ;
- aires au-dessus de la baseline ;
- ratio qualifiant / quantifiant.

Les composés partageant un même m/z sont départagés selon l'ordre fourni
par metadata.xlsx.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks, savgol_filter


def load_sample(sample_dir: str | Path):
    traces = {}
    for f in sorted(Path(sample_dir).glob("mz_*.csv")):
        m = re.fullmatch(r"mz_(\d+)(?:\.\d+)?", f.stem)
        if not m:
            continue
        df = pd.read_csv(f)
        traces[int(m.group(1))] = (
            df["time_min"].to_numpy(float),
            df["intensity"].to_numpy(float),
        )
    if not traces:
        raise FileNotFoundError(f"Aucun mz_*.csv dans {sample_dir}")
    return traces


def noise_sigma(y: np.ndarray) -> float:
    d2 = np.diff(y, n=2)
    mad = np.median(np.abs(d2 - np.median(d2)))
    return max(1.4826 * mad / np.sqrt(6.0), 1e-12)


def smooth(y: np.ndarray, window: int = 7) -> np.ndarray:
    window = min(window, len(y) - (1 - len(y) % 2))
    if window < 5:
        return y.copy()
    return savgol_filter(y, window | 1, 2)


def candidate_peaks(ys: np.ndarray, sigma: float, min_snr: float = 5.0):
    idx, props = find_peaks(ys, prominence=min_snr * sigma)
    return idx, props["prominences"]


def _descend(ys, start, step, sigma, lookahead, flat_k, rise_k, max_pts):
    n = len(ys)
    i, i_min = start, start
    for _ in range(max_pts):
        j = i + step
        if j < 0 or j >= n:
            return i, "limite"
        i = j
        if ys[i] < ys[i_min]:
            i_min = i
        elif ys[i] - ys[i_min] > rise_k * sigma:
            return i_min, "vallee"

        k = i + step * lookahead
        if 0 <= k < n and ys[i] - ys[k] < flat_k * sigma and ys[i] - ys[i_min] <= flat_k * sigma:
            far = min(max(i + step * 3 * lookahead, 0), n - 1)
            ahead = ys[min(i, far): max(i, far) + 1]
            if ahead.max() - ys[i] <= rise_k * sigma:
                return i, "baseline"
    return i, "limite"


@dataclass
class Peak:
    rt: float
    left: float
    right: float
    area: float
    snr: float
    left_stop: str
    right_stop: str


def integrate_peak(t, y, apex, *, ys=None, sigma=None, lookahead_min=0.03,
                   flat_k=1.0, rise_k=3.0, max_halfwidth_min=1.0):
    ys = smooth(y) if ys is None else ys
    sigma = noise_sigma(y) if sigma is None else sigma
    dt = float(np.median(np.diff(t)))
    look = max(3, int(round(lookahead_min / dt)))
    max_pts = max(look + 1, int(round(max_halfwidth_min / dt)))

    yd = uniform_filter1d(y, 5, mode="nearest")
    lo, hi = max(apex - 2, 0), min(apex + 2, len(yd) - 1)
    top = lo + int(np.argmax(yd[lo:hi + 1]))

    il, stop_l = _descend(yd, top, -1, sigma, look, flat_k, rise_k, max_pts)
    ir, stop_r = _descend(yd, top, +1, sigma, look, flat_k, rise_k, max_pts)
    il, ir = min(il, apex), max(ir, apex)

    level = lambda i: float(np.median(y[max(0, i - 2): i + 3]))
    baseline = min(level(il), level(ir))
    seg_t, seg_y = t[il:ir + 1], y[il:ir + 1]
    area = float((getattr(np, "trapezoid", None) or np.trapz)(seg_y - baseline, seg_t))

    rt = float(t[apex])
    if 0 < apex < len(ys) - 1:
        a, b, c = ys[apex - 1], ys[apex], ys[apex + 1]
        den = a - 2 * b + c
        if den < 0:
            rt += 0.5 * (a - c) / den * dt

    snr = float((ys[apex] - baseline) / sigma)
    return Peak(rt, float(t[il]), float(t[ir]), area, snr, stop_l, stop_r)


def compute_reference(traces, compounds: pd.DataFrame, min_snr: float = 5.0) -> pd.DataFrame:
    """compounds : name, mz_quant, mz_qual, elution_rank."""
    rows = []
    cache = {}

    def prep(mz):
        if mz not in cache:
            t, y = traces[mz]
            cache[mz] = (t, y, smooth(y), noise_sigma(y))
        return cache[mz]

    for mzq, grp in compounds.groupby("mz_quant", sort=False):
        grp = grp.sort_values("elution_rank")
        if int(mzq) not in traces:
            for _, c in grp.iterrows():
                rows.append({"name": c["name"], "warning": f"m/z quantifiant {mzq} absent"})
            continue
        t, y, ys, sig = prep(int(mzq))
        idx, prom = candidate_peaks(ys, sig, min_snr)

        n = len(grp)
        if len(idx) < n:
            for _, c in grp.iterrows():
                rows.append({"name": c["name"], "warning": f"{len(idx)} pic(s) pour {n} composé(s)"})
            continue

        order = np.argsort(prom)[::-1]
        keep = np.sort(idx[order[:n]])

        for apex, (_, c) in zip(keep, grp.iterrows()):
            pq = integrate_peak(t, y, apex, ys=ys, sigma=sig)

            if pd.isna(c["mz_qual"]) or int(c["mz_qual"]) not in traces:
                rows.append({"name": c["name"], "rt_ref": pq.rt, "area_quant": pq.area,
                             "warning": "m/z qualifiant absent ; ratio non calculable"})
                continue
            tl, yl, yls, sigl = prep(int(c["mz_qual"]))
            lo, hi = np.searchsorted(tl, [pq.left, pq.right])
            hi = min(hi, len(tl) - 1)
            if lo >= len(tl) or hi <= lo:
                rows.append({"name": c["name"], "rt_ref": pq.rt, "area_quant": pq.area,
                             "warning": "fenêtre qualifiante absente ou trop courte"})
                continue
            apex_l = lo + int(np.argmax(yls[lo:hi + 1]))
            pl = integrate_peak(tl, yl, apex_l, ys=yls, sigma=sigl)

            warning = []
            if "limite" in (pq.left_stop, pq.right_stop):
                warning.append("borne quantifiant à vérifier")
            if "limite" in (pl.left_stop, pl.right_stop):
                warning.append("borne qualifiant à vérifier")

            rows.append({
                "name": c["name"],
                "mz_quant": int(mzq),
                "mz_qual": int(c["mz_qual"]),
                "rt_ref": pq.rt,
                "quant_left": pq.left,
                "quant_right": pq.right,
                "qual_left": pl.left,
                "qual_right": pl.right,
                "area_quant": pq.area,
                "area_qual": pl.area,
                "ratio_ref": pl.area / pq.area if pq.area > 0 else np.nan,
                "snr_quant": pq.snr,
                "snr_qual": pl.snr,
                "warning": " ; ".join(warning),
            })

    if not rows:
        return pd.DataFrame(columns=["name", "elution_rank", "rt_ref", "ratio_ref", "warning"])
    return pd.DataFrame(rows).merge(
        compounds[["name", "elution_rank"]], on="name"
    ).sort_values("elution_rank").reset_index(drop=True)

