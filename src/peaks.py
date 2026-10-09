"""Traitement d'un pic chromatographique : bruit, apex, bornes et aire.

Un chromatogramme par ion est une série temporelle (time_min, intensity) : du
bruit autour d'une ligne de base, et des pics en cloche, souvent avec une
traînée à droite (les HAP « collent » à la paroi de la colonne).

Démarche, inspirée du travail manuel du chimiste :
  1. estimer le bruit de la trace, pour exprimer tous les seuils en multiples
     du bruit plutôt qu'en valeurs absolues ;
  2. repérer le sommet (apex) : c'est le temps de rétention ;
  3. descendre à gauche et à droite du sommet tant que le signal décroît
     (observation de la dérivée) : les points d'arrêt sont les bornes ;
  4. intégrer le signal brut entre les bornes, au-dessus de la ligne de base
     (somme de trapèzes, l'équivalent informatique de l'intégrale).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks, savgol_filter

_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy 1.x et 2.x


def noise_sigma(y: np.ndarray) -> float:
    """Écart-type du bruit, estimé sur les différences secondes (insensible aux pics
    et à la dérive lente). Var(diff2 d'un bruit blanc) = 6 sigma²."""
    d2 = np.diff(y, n=2)
    mad = np.median(np.abs(d2 - np.median(d2)))
    return max(1.4826 * mad / np.sqrt(6.0), 1e-12)


def smooth(y: np.ndarray, window: int = 7) -> np.ndarray:
    """Lissage Savitzky-Golay : sert à localiser (apex, bornes), jamais à intégrer."""
    window = min(window, len(y) - (1 - len(y) % 2))
    if window < 5:
        return y.copy()
    return savgol_filter(y, window | 1, 2)


def candidate_peaks(ys: np.ndarray, sigma: float, min_snr: float = 10.0):
    """Indices des maxima dont la proéminence dépasse min_snr x bruit."""
    idx, props = find_peaks(ys, prominence=min_snr * sigma)
    return idx, props["prominences"]


def _descend(ys, start, step, sigma, lookahead, flat_k, rise_k, max_pts):
    """Descend depuis `start` dans le sens `step` (+1 droite, -1 gauche).

    Arrêt sur le premier de ces critères :
      - 'vallee'   : le signal remonte de plus de rise_k x bruit au-dessus du
                     minimum courant -> pic voisin (coélution), borne = le minimum ;
      - 'baseline' : sur `lookahead` points le signal ne descend plus de plus de
                     flat_k x bruit -> retour à la ligne de base ;
      - 'limite'   : bord du signal ou largeur maximale atteinte.
    """
    n = len(ys)
    i, i_min = start, start
    for _ in range(max_pts):
        j = i + step
        if j < 0 or j >= n:
            return i, "fin_signal"       # bord de la fenêtre d'acquisition SIM
        i = j
        if ys[i] < ys[i_min]:
            i_min = i
        elif ys[i] - ys[i_min] > rise_k * sigma:
            return i_min, "vallee"
        k = i + step * lookahead
        if 0 <= k < n and ys[i] - ys[k] < flat_k * sigma and ys[i] - ys[i_min] <= flat_k * sigma:
            # Plat : vraie ligne de base, ou fond de vallée avant un pic voisin ?
            far = min(max(i + step * 3 * lookahead, 0), n - 1)
            ahead = ys[min(i, far): max(i, far) + 1]
            if ahead.max() - ys[i] <= rise_k * sigma:
                return i, "baseline"
    return i, "limite"


@dataclass
class Peak:
    rt: float            # temps de rétention de l'apex (min)
    left: float          # borne gauche (min)
    right: float         # borne droite (min)
    area: float          # aire au-dessus de la baseline
    height: float        # hauteur de l'apex au-dessus de la baseline
    baseline: float      # niveau de baseline à l'apex
    left_stop: str       # 'baseline' | 'vallee' | 'limite'
    right_stop: str
    snr: float


def integrate_peak(t, y, apex, *, ys=None, sigma=None, baseline="min_bornes", bounds=None,
                   lookahead_min=0.03, flat_k=1.0, rise_k=3.0, max_halfwidth_min=1.0) -> Peak:
    """Bornes par double descente depuis l'apex, puis aire sur le signal brut.

    bounds = (t_gauche, t_droite) impose les bornes au lieu de les chercher
             (utilisé pour intégrer le qualifiant sur la fenêtre du quantifiant).

    baseline = 'min_bornes' : horizontale au plus bas des deux bornes (convention de
               l'application de démonstration) ; 'lineaire' : droite entre les bornes.
    """
    t, y = np.asarray(t, float), np.asarray(y, float)
    if t.ndim != 1 or y.ndim != 1 or len(t) != len(y) or len(t) < 3:
        raise ValueError("Une trace doit contenir au moins trois couples temps/intensité")
    if not np.isfinite(t).all() or not np.isfinite(y).all() or not (np.diff(t) > 0).all():
        raise ValueError("Temps strictement croissants et valeurs finies requis")
    if not 0 <= apex < len(t):
        raise ValueError("Apex hors de la trace")
    ys = smooth(y) if ys is None else ys
    sigma = noise_sigma(y) if sigma is None else sigma
    dt = float(np.median(np.diff(t)))
    look = max(3, int(round(lookahead_min / dt)))
    max_pts = max(look + 1, int(round(max_halfwidth_min / dt)))

    if bounds is None:
        # Moyenne glissante pour la descente : contrairement à Savitzky-Golay, elle ne
        # crée pas de creux artificiel au pied du front montant.
        yd = uniform_filter1d(y, 5, mode="nearest")
        top = apex                      # sommet du signal de descente, voisin de l'apex
        # (±2 points au plus : au-delà, on grimperait sur un pic voisin plus haut)
        lo, hi = max(apex - 2, 0), min(apex + 2, len(yd) - 1)
        top = lo + int(np.argmax(yd[lo:hi + 1]))
        il, stop_l = _descend(yd, top, -1, sigma, look, flat_k, rise_k, max_pts)
        ir, stop_r = _descend(yd, top, +1, sigma, look, flat_k, rise_k, max_pts)
        il, ir = min(il, apex), max(ir, apex)
    else:
        il, ir = np.searchsorted(t, bounds[0]), min(np.searchsorted(t, bounds[1]), len(t) - 1)
        if il >= ir or not il <= apex <= ir:
            raise ValueError("Bornes sans intervalle valide autour de l’apex")
        stop_l = stop_r = "impose"

    # Niveau aux bornes : médiane locale du brut, moins sensible à un point bruité.
    lvl = lambda i: float(np.median(y[max(0, i - 2): i + 3]))
    yl, yr = lvl(il), lvl(ir)
    seg_t, seg_y = t[il: ir + 1], y[il: ir + 1]
    if baseline == "lineaire":
        base = np.interp(seg_t, [t[il], t[ir]], [yl, yr])
    elif baseline == "min_bornes":
        base = np.full_like(seg_t, min(yl, yr))
    else:
        raise ValueError(f"baseline inconnue : {baseline}")
    area = float(_trapz(seg_y - base, seg_t))

    # Apex affiné par parabole sur 3 points du signal lissé (sous-échantillon).
    rt = float(t[apex])
    if 0 < apex < len(ys) - 1:
        a, b, c = ys[apex - 1], ys[apex], ys[apex + 1]
        den = a - 2 * b + c
        if den < 0:
            rt += 0.5 * (a - c) / den * dt
    b_apex = float(base[apex - il])
    height = float(ys[apex] - b_apex)
    return Peak(rt, float(t[il]), float(t[ir]), area, height, b_apex,
                stop_l, stop_r, height / sigma)

