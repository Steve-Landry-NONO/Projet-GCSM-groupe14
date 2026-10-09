"""Application de la référence GAM-6 aux autres échantillons (GAM, SF).

Pour chaque échantillon, le chimiste cherche chaque composé près de son temps
de référence, vérifie son identité (TR ±0,2 %, ratio ±23 %) et calcule sa
réponse normalisée par l'étalon interne : c'est cette réponse qui alimente la
calibration.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from alerts import EDGE_MAX_RATIO, VALLEY_MAX_RATIO, Alert, blocking_from_reference, render
from peaks import candidate_peaks
from reference import RATIO_TOL_REL, RT_TOL_REL, _preparer, measure_peak


def measure_sample(traces, ref: pd.DataFrame, compounds: pd.DataFrame, *, search_rel=0.01,
                   min_snr=3.0, sample_type=None, baseline="min_bornes", qual_bounds="quantifiant",
                   valley_max=VALLEY_MAX_RATIO, edge_max=EDGE_MAX_RATIO, **kw) -> pd.DataFrame:
    """Mesure chaque composé d'un échantillon à partir de la référence GAM-6.

    Règle métier : on prend le pic du m/z quantifiant le plus proche du temps de
    référence (recherche dans ±search_rel), puis on contrôle
      - le temps de rétention : |TR - TR_ref| <= 0,2 % de TR_ref ;
      - le ratio qualifiant/quantifiant : |ratio - ratio_ref| <= 23 % de ratio_ref ;
    et on calcule la réponse normalisée aire_quant / aire_quant(ISTD).

    Une mesure est « OK » si TR et ratio passent et qu'aucune alerte bloquante ne
    subsiste (voir alerts.py) ; les alertes informatives restent affichées.
    """
    prep = _preparer(traces)
    istd_of = dict(zip(compounds["name"], compounds["istd"].fillna("")))
    rows = []
    expected = compounds.merge(ref, on="name", how="left", suffixes=("_method", ""))
    for _, r in expected.iterrows():
        nm = r["name"]
        mzq = int(r.get("mz_quant_method", r.get("mz_quant")))
        mzl = r.get("mz_qual_method", r.get("mz_qual"))
        rt0 = r.get("rt_ref", np.nan)
        row = {"name": nm, "istd": istd_of.get(nm, ""), "rt_ref": rt0, "ratio_ref": r.get("ratio_ref")}
        use = str(r.get("use", "")).strip().lower()
        applicable = str(r.get("applicable_to", ""))
        if use == "inutilisé" or (sample_type and applicable not in ("", "nan")
                                  and sample_type not in {v.strip() for v in applicable.split(",")}):
            rows.append({**row, "status": "non applicable", "warning": "exclu par les métadonnées"})
            continue
        if not np.isfinite(rt0):
            rows.append({**row, "status": "référence absente", "warning": "temps de référence manquant"})
            continue
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
                                    baseline=baseline, qual_bounds=qual_bounds,
                                    valley_max=valley_max, edge_max=edge_max, **kw)
        d_rt = (pq.rt - rt0) / rt0
        row.update(rt=pq.rt, d_rt_pct=100 * d_rt, rt_ok=bool(abs(d_rt) <= RT_TOL_REL),
                   left=pq.left, right=pq.right, area_quant=pq.area, snr_quant=pq.snr,
                   valley_ratio=pq.bound_ratio("vallee"), edge_ratio=pq.bound_ratio("fin_signal"))
        if pl is not None and pq.area > 0:
            ratio = pl.area / pq.area
            ratio0 = r.get("ratio_ref", np.nan)
            d_ratio = (ratio - ratio0) / ratio0 if pd.notna(ratio0) and ratio0 > 0 else np.nan
            row.update(area_qual=pl.area, ratio=ratio, d_ratio_pct=100 * d_ratio,
                       ratio_ok=bool(abs(d_ratio) <= RATIO_TOL_REL))
        warn += blocking_from_reference(r)
        if r.get("double_peak") is True or r.get("double_peak") == 1:
            warn.append(Alert("double pic (double_peak) : règle métier non définie, aire non validée"))
        if pd.notna(r.get("median_despike")) and r.get("median_despike") == 1:
            warn.append(Alert("filtrage des pics parasites (median_despike) non implémenté"))
        if not row["rt_ok"]:
            warn.append(Alert("TR hors ±0,2 %"))
        if row.get("ratio_ok") is not True:
            warn.append(Alert("ratio hors ±23 %" if row.get("ratio_ok") is False
                              else "ratio qualifiant/quantifiant non calculable"))
        row.update(render(warn))
        row["status"] = "OK" if not row["blocking_alerts"] else "à vérifier"
        rows.append(row)
    out = pd.DataFrame(rows)
    # On conserve toutes les lignes, même lorsque rien n'a été détecté.
    for col in ("area_quant", "rt", "ratio", "rt_ok", "ratio_ok", "warning", "blocking_alerts", "info_alerts"):
        if col not in out:
            out[col] = "" if col in ("warning", "blocking_alerts", "info_alerts") else np.nan
    for col in ("warning", "blocking_alerts", "info_alerts"):
        out[col] = out[col].fillna("")
    out["valley_max"], out["edge_max"] = valley_max, edge_max

    def add_blocking(mask, message):
        for col in ("warning", "blocking_alerts"):
            out.loc[mask, col] = (out.loc[mask, col] + " ; " + message).str.strip(" ;")
        out.loc[mask & out["status"].eq("OK"), "status"] = "à vérifier"

    # Deux cibles ne doivent pas utiliser le même sommet sur un ion partagé.
    mz_of = dict(zip(compounds["name"], compounds["mz_quant"]))
    out["mz_quant"] = out["name"].map(mz_of)
    duplicated = out["rt"].notna() & out.duplicated(["mz_quant", "rt"], keep=False)
    add_blocking(duplicated, "apex attribué à plusieurs composés")
    area = out.set_index("name")["area_quant"].to_dict()
    status = out.set_index("name")["status"].to_dict()
    out["area_istd"] = out["istd"].map(area)
    out["istd_status"] = out["istd"].map(status)
    valid = (out["status"].eq("OK") & out["istd_status"].eq("OK")
             & np.isfinite(out["area_quant"]) & np.isfinite(out["area_istd"])
             & out["area_quant"].gt(0) & out["area_istd"].gt(0))
    out["response"] = np.nan
    out.loc[valid, "response"] = out.loc[valid, "area_quant"] / out.loc[valid, "area_istd"]
    # Une réponse exploratoire conserve le calcul sans certifier son aire.
    out["identity_ok"] = out["rt_ok"].eq(True) & out["ratio_ok"].eq(True) & ~duplicated
    identity = out.set_index("name")["identity_ok"].to_dict()
    out["istd_identity_ok"] = out["istd"].map(identity).eq(True)
    exploratory = (out["identity_ok"] & out["istd_identity_ok"]
                   & np.isfinite(out["area_quant"]) & np.isfinite(out["area_istd"])
                   & out["area_quant"].gt(0) & out["area_istd"].gt(0))
    out["response_exploratory"] = np.nan
    out.loc[exploratory, "response_exploratory"] = out.loc[exploratory, "area_quant"] / out.loc[exploratory, "area_istd"]
    has_istd = out["istd"].ne("")
    istd_bad = has_istd & ~out["istd_status"].eq("OK")
    for col in ("warning", "blocking_alerts"):
        out.loc[istd_bad, col] = (out.loc[istd_bad, col] + " ; étalon interne non validé").str.strip(" ;")
    out.loc[istd_bad & out["status"].eq("OK"), "status"] = "ISTD à vérifier"
    return out
