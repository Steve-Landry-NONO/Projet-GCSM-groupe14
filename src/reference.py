"""Référence chromatographique (GAM-6) : temps de rétention et ratios de référence.

Pour chaque composé de la méthode (metadata.xlsx) :
  1. repérer son pic sur le m/z quantifiant ; si plusieurs composés partagent
     ce m/z, les départager par l'ordre d'élution, repère de la méthode, à confirmer ;
  2. apex -> temps de rétention de référence ;
  3. bornes par descente à gauche et à droite du sommet ;
  4. aire au-dessus de la ligne de base, sur le signal brut ;
  5. même traitement sur le m/z qualifiant, ancré sur le pic quantifiant ;
  6. ratio de référence = aire qualifiant / aire quantifiant.

Les seuils sont exprimés en multiples du bruit de chaque trace : aucun temps de
rétention absolu n'est codé en dur, car il dépend de l'état de la colonne.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from peaks import candidate_peaks, integrate_peak, noise_sigma, smooth

# Tolérances des notes : paramètres de travail à confirmer avec le professeur.
RT_TOL_REL = 0.002     # ±0,2 % sur le temps de rétention
RATIO_TOL_REL = 0.23   # ±23 % sur le ratio qualifiant/quantifiant


def load_sample(sample_dir: str | Path) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """Lit tous les mz_XXX.csv d'un échantillon -> {m/z: (time_min, intensity)}."""
    traces = {}
    for f in sorted(Path(sample_dir).glob("mz_*.csv")):
        m = re.fullmatch(r"mz_(\d+)(?:\.\d+)?", f.stem)
        if not m:
            continue
        df = pd.read_csv(f)
        traces[int(m.group(1))] = (df["time_min"].to_numpy(float),
                                   df["intensity"].to_numpy(float))
    if not traces:
        raise FileNotFoundError(f"Aucun mz_*.csv dans {sample_dir}")
    return traces


def default_method(root: Path) -> Path:
    """Emplacement par défaut de metadata.xlsx : data/ puis la racine du projet."""
    for cand in (root / "data" / "metadata.xlsx", root / "metadata.xlsx"):
        if cand.exists():
            return cand
    return root / "data" / "metadata.xlsx"


def load_compounds(path: str | Path) -> pd.DataFrame:
    """Table des composés depuis metadata.xlsx (source de vérité) ou un CSV équivalent.

    Dans metadata.xlsx, l'ordre des lignes de `compound_informations` est l'ordre
    d'élution ; `compound_deuterated` associe chaque cible à son étalon interne.
    """
    path = Path(path)
    if path.suffix.lower() in (".xlsx", ".xls"):
        info = pd.read_excel(path, sheet_name="compound_informations")
        deut = pd.read_excel(path, sheet_name="compound_deuterated")
        c = info.rename(columns={"compound_name": "name", "quantifier": "mz_quant",
                                 "qualifier": "mz_qual"})
        c = c.dropna(subset=["name", "mz_quant"]).reset_index(drop=True)
        c["elution_rank"] = np.arange(1, len(c) + 1)
        c["istd"] = c["name"].map(dict(zip(deut["compound_name"], deut["deuterated_name"])))
    else:
        c = pd.read_csv(path)
    c["mz_quant"] = c["mz_quant"].astype(int)
    c["mz_qual"] = pd.to_numeric(c["mz_qual"], errors="coerce")
    if "istd" not in c:
        c["istd"] = np.nan
    return c


def measure_peak(prep, traces, apex, mzq, mzl, *, min_snr=10.0, baseline="min_bornes",
                 qual_bounds="quantifiant", **kw):
    """Intègre le pic quantifiant à `apex`, puis le qualifiant ancré dessus.

    Retourne (pic quantifiant, pic qualifiant ou None, liste d'alertes).
    """
    warn = []
    t, y, ys, sig = prep(mzq)
    pq = integrate_peak(t, y, apex, ys=ys, sigma=sig, baseline=baseline, **kw)
    if "limite" in (pq.left_stop, pq.right_stop):
        warn.append("borne non trouvée (largeur maximale atteinte)")
    if "fin_signal" in (pq.left_stop, pq.right_stop):
        warn.append("traînée coupée par la fin de la fenêtre SIM")
    if pd.isna(mzl) or int(mzl) not in traces:
        warn.append(f"m/z qualifiant {mzl} absent")
        return pq, None, warn
    tl, yl, yls, sigl = prep(int(mzl))
    # Apex qualifiant : maximum local au plus près du temps du quantifiant
    # (le m/z qualifiant peut porter le pic d'un autre composé dans la fenêtre).
    dtl = float(np.median(np.diff(tl)))
    half = max(3 * dtl, pq.rt * RT_TOL_REL)
    lo, hi = np.searchsorted(tl, [pq.rt - half, pq.rt + half])
    hi = min(hi, len(tl) - 1)
    if lo >= len(tl) or hi <= lo:
        warn.append("fenêtre qualifiante absente ou trop courte")
        return pq, None, warn
    apex_l = lo + int(np.argmax(yls[lo:hi + 1]))
    own = integrate_peak(tl, yl, apex_l, ys=yls, sigma=sigl, baseline=baseline, **kw)
    if qual_bounds == "quantifiant":
        # Bornes du quantifiant, resserrées aux vallées du qualifiant s'il y a
        # un pic voisin sur ce m/z.
        left = max(pq.left, own.left) if own.left_stop == "vallee" else pq.left
        right = min(pq.right, own.right) if own.right_stop == "vallee" else pq.right
        pl = integrate_peak(tl, yl, apex_l, ys=yls, sigma=sigl, baseline=baseline,
                            bounds=(left, right), **kw)
    else:
        pl = own
    if abs(pl.rt - pq.rt) > pq.rt * RT_TOL_REL:
        warn.append("apex qualifiant décalé de plus de 0,2 %")
    if pl.snr < min_snr:
        warn.append("qualifiant faible (S/B < seuil)")
    return pq, pl, warn


def _preparer(traces):
    cache = {}

    def prep(mz):
        if mz not in cache:
            t, y = traces[mz]
            cache[mz] = (t, y, smooth(y), noise_sigma(y))
        return cache[mz]
    return prep


def compute_reference(traces, compounds: pd.DataFrame, *, min_snr=10.0,
                      baseline="min_bornes", qual_bounds="quantifiant", **kw) -> pd.DataFrame:
    """Temps de rétention et ratio de référence pour chaque composé.

    compounds : colonnes `name`, `mz_quant`, `mz_qual`, `elution_rank`
        (`elution_rank` = ordre de sortie de colonne, repère de la méthode, à confirmer ; il
        départage les composés qui partagent un même m/z quantifiant).
    qual_bounds : 'quantifiant' = le qualifiant est intégré sur les bornes du
        quantifiant ; 'independant' = il a sa propre descente.
    """
    prep = _preparer(traces)

    comps = compounds.sort_values("elution_rank").reset_index(drop=True)
    picks, warns = {}, {n: [] for n in comps["name"]}

    # 1) Affectation par m/z : les n pics les plus nets, rangés dans l'ordre d'élution.
    for mzq, grp in comps.groupby("mz_quant", sort=False):
        if mzq not in traces:
            for n in grp["name"]:
                warns[n].append(f"m/z {mzq} absent")
            continue
        t, y, ys, sig = prep(mzq)
        idx, prom = candidate_peaks(ys, sig, min_snr)
        n = len(grp)
        if len(idx) < n:
            for nm in grp["name"]:
                warns[nm].append(f"{len(idx)} pic(s) pour {n} composé(s) sur m/z {mzq}")
            continue
        order = np.argsort(prom)[::-1]
        keep = np.sort(idx[order[:n]])
        ambiguous = len(idx) > n and prom[order[n]] > 0.2 * prom[order[n - 1]]
        for apex, nm in zip(keep, grp["name"]):
            picks[nm] = int(apex)
            if ambiguous:
                warns[nm].append("autre pic proche en intensité sur ce m/z")

    # 2) Contrôle croisé par l'ordre d'élution : un composé doit sortir entre ses
    #    voisins (rangs inférieur et supérieur) trouvés sur les autres m/z. Sinon on
    #    reprend le pic le plus net de son m/z dans cette fenêtre. Cas typique : un
    #    deutéré dont le m/z est partagé avec un surrogat absent de la table.
    rt_of = lambda nm: prep(int(comps.loc[comps["name"] == nm, "mz_quant"].iloc[0]))[0][picks[nm]]
    for i, c in comps.iterrows():
        nm = c["name"]
        if nm not in picks:
            continue
        before = [rt_of(o) for o in comps["name"][:i] if o in picks]
        after = [rt_of(o) for o in comps["name"][i + 1:] if o in picks]
        lo, hi = (max(before) if before else -np.inf), (min(after) if after else np.inf)
        if lo < rt_of(nm) < hi or lo >= hi:
            continue
        t, y, ys, sig = prep(int(c["mz_quant"]))
        idx, prom = candidate_peaks(ys, sig, min_snr)
        ok = (t[idx] > lo) & (t[idx] < hi)
        if ok.any():
            picks[nm] = int(idx[ok][np.argmax(prom[ok])])
            warns[nm].append("pic réaffecté d'après l'ordre d'élution")

    # 3) Intégration quantifiant + qualifiant.
    rows = []
    for _, c in comps.iterrows():
        nm, mzq, mzl = c["name"], int(c["mz_quant"]), c["mz_qual"]
        row = {"elution_rank": c["elution_rank"], "name": nm, "mz_quant": mzq, "mz_qual": mzl}
        warn = warns[nm]
        if nm not in picks:
            rows.append({**row, "warning": " ; ".join(warn)})
            continue
        pq, pl, w = measure_peak(prep, traces, picks[nm], mzq, mzl, min_snr=min_snr,
                                 baseline=baseline, qual_bounds=qual_bounds, **kw)
        warn += w
        if "vallee" in (pq.left_stop, pq.right_stop):
            warn.append("pic voisin : séparation des aires à valider")
        row.update(rt_ref=pq.rt, rt_min=pq.rt * (1 - RT_TOL_REL), rt_max=pq.rt * (1 + RT_TOL_REL),
                   q_left=pq.left, q_right=pq.right, q_stop=f"{pq.left_stop}/{pq.right_stop}",
                   area_quant=pq.area, snr_quant=pq.snr)
        if pl is not None:
            ratio = pl.area / pq.area if pq.area > 0 else np.nan
            row.update(l_left=pl.left, l_right=pl.right, area_qual=pl.area, snr_qual=pl.snr,
                       d_rt_qual=pl.rt - pq.rt, ratio_ref=ratio,
                       ratio_min=ratio * (1 - RATIO_TOL_REL), ratio_max=ratio * (1 + RATIO_TOL_REL))
        row["warning"] = " ; ".join(warn)
        rows.append(row)

    out = pd.DataFrame(rows)
    for col in ("rt_ref", "ratio_ref", "area_quant", "area_qual"):
        if col not in out:
            out[col] = np.nan
    # Garde-fou final : les temps de référence doivent croître avec l'ordre d'élution.
    if "rt_ref" in out:
        rt = out["rt_ref"].to_numpy(float)
        prev = np.fmax.accumulate(np.r_[-np.inf, rt[:-1]])
        bad = rt <= prev
        out.loc[bad, "warning"] = (out.loc[bad, "warning"].fillna("")
                                   + " ; ordre d'élution non respecté").str.strip(" ;")
    return out

