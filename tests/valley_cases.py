"""Cas synthétiques de deux pics voisins sur un même ion.

Chaque cas combine une proportion de hauteurs, un écart entre apex, un niveau de
bruit et une ligne de base. La vérité terrain (aire de chaque pic) est connue
par construction, et un témoin « pic seul » sépare l'erreur due au voisin de
celle due au bruit.

Utilisé par tests/test_alerts.py (assertions) et par
tests/caracteriser_vallee.py (tableau de caractérisation documenté).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.stats import exponnorm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from peaks import candidate_peaks, integrate_peak, noise_sigma, smooth  # noqa: E402

DT, SIGMA, TAU = 0.0037, 0.008, 0.03       # pas réel du batch, largeur et traînée
T = np.arange(10, 12, DT)
BASELINES = {
    "plate": lambda h: np.zeros_like(T),
    "pente": lambda h: 0.02 * h * (T - T[0]),
    "decalee": lambda h: 0.05 * h + 0.01 * h * (T - T[0]),
}


def emg(loc: float, height: float) -> np.ndarray:
    """Pic gaussien à traînée exponentielle, de hauteur maximale `height`."""
    p = exponnorm.pdf(T, TAU / SIGMA, loc=loc, scale=SIGMA)
    return height * p / p.max()


def measure(y: np.ndarray, n: int):
    """Les n pics les plus nets, intégrés comme dans le pipeline."""
    ys, sg = smooth(y), noise_sigma(y)
    idx, prom = candidate_peaks(ys, sg, 10)
    if len(idx) < n:
        return None
    top = np.sort(idx[np.argsort(prom)[::-1][:n]])
    return [integrate_peak(T, y, int(a), ys=ys, sigma=sg) for a in top]


def run_case(height_ratio, sep, noise, base, seed=0, with_peaks=False):
    """Mesure un couple de pics ; une ligne par pic.

    height_ratio : hauteur du pic 2 / hauteur du pic 1 ;
    sep : écart entre les positions des deux pics (min) ;
    noise : écart-type du bruit en fraction de la plus grande hauteur.
    """
    rng = np.random.default_rng(seed)
    h1, h2 = 1e5, 1e5 * height_ratio
    hmax = max(h1, h2)
    p1, p2 = emg(11.0, h1), emg(11.0 + sep, h2)
    areas = (np.trapezoid(p1, T), np.trapezoid(p2, T)) if hasattr(np, "trapezoid") \
        else (np.trapz(p1, T), np.trapz(p2, T))
    ref = measure(p1 + p2, 2)                    # sans bruit, ligne de base nulle
    y = p1 + p2 + BASELINES[base](hmax) + rng.normal(0, noise * hmax, T.size)
    got = measure(y, 2)
    rows = []
    for k, (alone, area) in enumerate(((p1, areas[0]), (p2, areas[1]))):
        iso = measure(alone + BASELINES[base](hmax) + rng.normal(0, noise * hmax, T.size), 1)
        row = dict(height_ratio=height_ratio, sep=sep, noise=noise, base=base, seed=seed,
                   peak=k + 1, detected=got is not None and ref is not None,
                   v_ref=ref[k].bound_ratio("vallee") if ref else np.nan,
                   err_iso=100 * (iso[0].area - area) / area if iso else np.nan)
        if got is not None:
            if with_peaks:
                row["peak_obj"] = got[k]
            row.update(v_est=got[k].bound_ratio("vallee"),
                       err_area=100 * (got[k].area - area) / area,
                       stops=f"{got[k].left_stop}/{got[k].right_stop}")
        rows.append(row)
    return rows


GRID = dict(height_ratio=(1.0, 0.2, 5.0, 0.05),
            sep=(0.045, 0.06, 0.08, 0.10, 0.12, 0.14, 0.18),
            noise=(0.001, 0.01, 0.03),
            base=tuple(BASELINES))
