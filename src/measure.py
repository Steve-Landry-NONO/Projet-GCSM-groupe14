"""Application de la référence GAM-6 aux autres échantillons (GAM, SF).

Pour chaque échantillon, le chimiste cherche chaque composé près de son temps
de référence, vérifie son identité (TR ±0,2 %, ratio ±23 %) et calcule sa
réponse normalisée par l'étalon interne : c'est cette réponse qui alimente la
calibration.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from peaks import candidate_peaks
from reference import RATIO_TOL_REL, RT_TOL_REL, _preparer, measure_peak


def measure_sample(traces, ref: pd.DataFrame, compounds: pd.DataFrame, *, search_rel=0.01,
                   min_snr=3.0, baseline="min_bornes", qual_bounds="quantifiant", **kw) -> pd.DataFrame:
    """Mesure chaque composé d'un échantillon à partir de la référence GAM-6.

    Règle métier : on prend le pic du m/z quantifiant le plus proche du temps de
    référence (recherche dans ±search_rel), puis on contrôle
      - le temps de rétention : |TR - TR_ref| <= 0,2 % de TR_ref ;
      - le ratio qualifiant/quantifiant : |ratio - ratio_ref| <= 23 % de ratio_ref ;
    et on calcule la réponse normalisée aire_quant / aire_quant(ISTD).
    """
    prep = _preparer(traces)
    istd_of = dict(zip(compounds["name"], compounds["istd"].fillna("")))
    rows = []
    for _, r in ref.dropna(subset=["rt_ref"]).iterrows():
        nm, mzq, mzl, rt0 = r["name"], int(r["mz_quant"]), r["mz_qual"], float(r["rt_ref"])
        row = {"name": nm, "istd": istd_of.get(nm, ""), "rt_ref": rt0, "ratio_ref": r.get("ratio_ref")}
        if mzq not in traces:
            rows.append({**row, "status": "m/z absent"})
            continue
        t, y, ys, sig = prep(mzq)
        idx, _ = candidate_peaks(ys, sig, min_snr)
        idx = idx[np.abs(t[idx] - rt0) <= rt0 * search_rel]
        if len(idx) == 0:
            rows.append({**row, "status": "non détecté"})
            continue
        apex = int(idx[np.argmin(np.abs(t[idx] - rt0))])
        pq, pl, warn = measure_peak(prep, traces, apex, mzq, mzl, min_snr=min_snr,
                                    baseline=baseline, qual_bounds=qual_bounds, **kw)
        d_rt = (pq.rt - rt0) / rt0
        row.update(rt=pq.rt, d_rt_pct=100 * d_rt, rt_ok=bool(abs(d_rt) <= RT_TOL_REL),
                   left=pq.left, right=pq.right, area_quant=pq.area, snr_quant=pq.snr)
        if pl is not None and pq.area > 0:
            ratio = pl.area / pq.area
            d_ratio = (ratio - r["ratio_ref"]) / r["ratio_ref"]
            row.update(area_qual=pl.area, ratio=ratio, d_ratio_pct=100 * d_ratio,
                       ratio_ok=bool(abs(d_ratio) <= RATIO_TOL_REL))
        checks = [row.get("rt_ok"), row.get("ratio_ok")]
        row["status"] = "OK" if all(c is True for c in checks) else "à vérifier"
        if not row["rt_ok"]:
            warn.append("TR hors ±0,2 %")
        if row.get("ratio_ok") is False:
            warn.append("ratio hors ±23 %")
        row["warning"] = " ; ".join(warn)
        rows.append(row)
    out = pd.DataFrame(rows)
    # Réponse normalisée par l'étalon interne (ce qui alimentera la calibration).
    area = dict(zip(out["name"], out.get("area_quant", pd.Series(dtype=float))))
    out["area_istd"] = out["istd"].map(lambda i: area.get(i, np.nan) if i else np.nan)
    out["response"] = out["area_quant"] / out["area_istd"]
    return out
