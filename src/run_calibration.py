"""Calibration quadratique par HAP (gamme GAM) et contrôle des SF.

Usage :
    python src/run_calibration.py                # lit outputs/<batch>/peaks.csv
    python src/run_calibration.py --sf-tol 0.20      # tolérance SF (provisoire, à valider)

Modèle (comme MassHunter) : réponse = a·C² + b·C + c, ajusté par moindres carrés
pondérés en 1/C sur les 8 points GAM, avec réponse = aire quantifiant / aire ISTD.
Inversion : racine positive du polynôme, dans le domaine de la gamme.

Sorties dans outputs/<batch>/ :
    calibrations.csv          a, b, c, R², écart max des points recalculés, domaine
    calibration_points.csv    chaque point GAM : réponse, concentration recalculée, écart
    sf_results.csv            chaque SF x composé : concentration calculée, écart, PASS/FAIL
    calibration_courbes.png   une courbe par HAP avec les points GAM et SF
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def fit_quadratic(conc, resp, weight="1/x"):
    """Moindres carrés pondérés : retourne (a, b, c)."""
    conc, resp = np.asarray(conc, float), np.asarray(resp, float)
    w = 1 / conc if weight == "1/x" else np.ones_like(conc)
    A = np.c_[conc**2, conc, np.ones_like(conc)]
    sw = np.sqrt(w)
    a, b, c = np.linalg.lstsq(A * sw[:, None], resp * sw, rcond=None)[0]
    return float(a), float(b), float(c)


def r_squared(conc, resp, coef):
    pred = np.polyval(coef, conc)
    total = np.sum((resp - np.mean(resp)) ** 2)
    return float(1 - np.sum((resp - pred) ** 2) / total) if total > 0 else np.nan


def invert(coef, response, c_max, c_min=0.0):
    """Inverse une courbe croissante ; toute extrapolation reste signalée."""
    if not np.isfinite(response):
        return np.nan, "réponse manquante"
    a, b, c = np.asarray(coef, float)
    if not np.isfinite([a, b, c, c_min, c_max]).all() or c_min < 0 or c_max <= c_min:
        return np.nan, "calibration invalide"
    if min(2 * a * c_min + b, 2 * a * c_max + b) <= 0:
        return np.nan, "courbe non croissante sur la gamme"
    roots = np.roots([a, b, c - response]) if abs(a) > 1e-12 else np.array([(response - c) / b])
    real = roots[np.abs(roots.imag) < 1e-9].real
    # La pente positive choisit la branche physique, pas la plus petite racine.
    real = real[(real >= 0) & (2 * a * real + b > 0)]
    if len(real) != 1:
        return np.nan, "solution absente ou ambiguë"
    conc = float(real[0])
    eps = 1e-9 * max(1.0, c_max)
    return conc, ("ok" if c_min - eps <= conc <= c_max + eps else "hors gamme")


def sf_verdict(peak_status, inversion, error_pct, tol_pct):
    """Un petit écart ne suffit pas si le pic ou la calibration est douteux."""
    if peak_status != "OK" or inversion != "ok" or not np.isfinite(error_pct):
        return "NON VALIDÉ"
    return "PASS" if abs(error_pct) <= tol_pct else "FAIL"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=ROOT / "outputs")
    ap.add_argument("--batch", default=None, help="sous-dossier de outputs (par défaut : le seul présent)")
    ap.add_argument("--sf-tol", type=float, default=0.20,
                    help="écart relatif accepté sur les SF (défaut 0,20 = ±20 %%, PROVISOIRE)")
    ap.add_argument("--weight", choices=["1/x", "aucune"], default="1/x")
    args = ap.parse_args()

    batches = [d for d in args.out.iterdir() if (d / "peaks.csv").exists()] if args.out.exists() else []
    if args.batch:
        batches = [args.out / args.batch]
    if len(batches) != 1:
        raise SystemExit(f"Précise --batch parmi : {[b.name for b in batches]}" if batches
                         else "Aucun peaks.csv : lance d'abord run_reference.py puis run_batch.py")
    out = batches[0]
    peaks = pd.read_csv(out / "peaks.csv")
    targets = peaks[peaks["istd"].notna() & (peaks["istd"] != "")]

    gam = targets[targets["type"] == "GAM"]
    cal_rows, pt_rows = [], []
    for name, d in gam.groupby("name", sort=False):
        all_points = d.copy()
        valid = (d["status"].eq("OK") & np.isfinite(d["response"])
                 & np.isfinite(d["conc_nominal_ppm"]) & d["conc_nominal_ppm"].gt(0))
        if "istd_status" in d:
            valid &= d["istd_status"].eq("OK")
        if "warning" in d:
            valid &= d["warning"].fillna("").str.strip().eq("")
        d = d[valid].sort_values("conc_nominal_ppm")
        for _, excluded in all_points.loc[~valid].iterrows():
            pt_rows.append({"name": name, "sample": excluded["sample"], "level": excluded["level"],
                            "conc_nominal_ppm": excluded["conc_nominal_ppm"], "response": excluded["response"],
                            "status_point": excluded["status"], "included": False,
                            "exclusion_reason": "mesure ou ISTD non validé"})
        x, y = d["conc_nominal_ppm"].to_numpy(float), d["response"].to_numpy(float)
        if d["conc_nominal_ppm"].nunique() < 4:
            for _, r in d.iterrows():
                pt_rows.append({"name": name, "sample": r["sample"], "level": r["level"],
                    "conc_nominal_ppm": r["conc_nominal_ppm"], "response": r["response"],
                    "status_point": r["status"], "included": False,
                    "exclusion_reason": "moins de quatre concentrations distinctes valides"})
            cal_rows.append({"name": name, "status": f"{d.conc_nominal_ppm.nunique()} concentrations valides seulement"})
            continue
        coef = fit_quadratic(x, y, args.weight)
        back = [invert(coef, r, x.max(), x.min())[0] for r in y]
        err = 100 * (np.array(back) - x) / x
        for (_, r), bc, e in zip(d.iterrows(), back, err):
            pt_rows.append({"name": name, "sample": r["sample"], "level": r["level"],
                            "conc_nominal_ppm": r["conc_nominal_ppm"], "response": r["response"],
                            "response_fit": float(np.polyval(coef, r["conc_nominal_ppm"])),
                            "conc_backcalc_ppm": bc, "ecart_pct": e,
                            "status_point": r["status"], "included": True})
        cal_rows.append({"name": name, "istd": d["istd"].iloc[0], "a": coef[0], "b": coef[1],
                         "c": coef[2], "r2": r_squared(x, y, coef), "n_points": len(d),
                         "ecart_max_pct": float(np.nanmax(np.abs(err))) if np.isfinite(err).any() else np.nan,
                         "c_min_ppm": x.min(), "c_max_ppm": x.max(), "ponderation": args.weight,
                         "status": "ok" if min(2 * coef[0] * x.min() + coef[1],
                                                   2 * coef[0] * x.max() + coef[1]) > 0
                                  else "courbe non croissante"})
    cal = pd.DataFrame(cal_rows).reindex(columns=["name", "istd", "a", "b", "c", "r2", "n_points",
        "ecart_max_pct", "c_min_ppm", "c_max_ppm", "ponderation", "status"])
    cal.to_csv(out / "calibrations.csv", index=False, float_format="%.6g")
    points = pd.DataFrame(pt_rows).reindex(columns=["name", "sample", "level", "conc_nominal_ppm",
        "response", "response_fit", "conc_backcalc_ppm", "ecart_pct", "status_point", "included", "exclusion_reason"])
    points.to_csv(out / "calibration_points.csv", index=False, float_format="%.6g")

    sf_rows = []
    coefs = {r["name"]: r for _, r in cal.iterrows() if r.get("status") == "ok"}
    for _, r in targets[targets["type"] == "SF"].iterrows():
        row = {"sample": r["sample"], "name": r["name"], "response": r["response"],
               "conc_nominal_ppm": r["conc_nominal_ppm"], "peak_status": r["status"]}
        k = coefs.get(r["name"])
        if k is None:
            sf_rows.append({**row, "sf_status": "pas de calibration"})
            continue
        conc, why = invert((k["a"], k["b"], k["c"]), r["response"], k["c_max_ppm"], k["c_min_ppm"])
        mult = float(r["multiplier"]) if "multiplier" in r and pd.notna(r["multiplier"]) else 1.0
        final = conc * mult
        err = (final - r["conc_nominal_ppm"]) / r["conc_nominal_ppm"]
        status = sf_verdict(r["status"], why, 100 * err, 100 * args.sf_tol)
        sf_rows.append({**row, "conc_calc_ppm": conc, "multiplier": mult, "conc_final_ppm": final,
                        "ecart_pct": 100 * err, "inversion": why, "sf_status": status})
    sf = pd.DataFrame(sf_rows).reindex(columns=["sample", "name", "response", "conc_nominal_ppm",
        "peak_status", "conc_calc_ppm", "multiplier", "conc_final_ppm", "ecart_pct", "inversion", "sf_status"])
    sf.to_csv(out / "sf_results.csv", index=False, float_format="%.6g")

    plot(cal, points, sf, out / "calibration_courbes.png")

    with pd.option_context("display.width", 200):
        print(cal[["name", "a", "b", "c", "r2", "ecart_max_pct"]].to_string(
            index=False, float_format=lambda v: f"{v:.4f}"))
        if len(sf):
            print(f"\nContrôle SF (tolérance ±{100 * args.sf_tol:.0f} %, à valider) :")
            piv = sf.pivot_table(index="name", columns="sample", values="ecart_pct")
            print(piv.round(1).to_string())
            print("\n", sf["sf_status"].value_counts().to_string())
    print(f"\nRésultats : {out / 'calibrations.csv'}, {out / 'sf_results.csv'}")


def plot(cal, pts, sf, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ok = cal[cal.get("status", "ok") == "ok"]
    cols = 4
    if ok.empty:
        path.unlink(missing_ok=True)
        return
    rows = -(-len(ok) // cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.0 * rows), squeeze=False)
    for ax, (_, k) in zip(axes.flat, ok.iterrows()):
        p = pts[(pts["name"] == k["name"]) & pts["included"].eq(True)]
        xs = np.linspace(0, k["c_max_ppm"], 200)
        ax.plot(xs, np.polyval([k["a"], k["b"], k["c"]], xs), color="black", lw=0.9)
        ax.plot(p["conc_nominal_ppm"], p["response"], "o", mfc="white", mec="black", ms=4, label="GAM")
        s = sf[sf["name"] == k["name"]] if len(sf) else sf
        if len(s):
            ax.plot(s["conc_calc_ppm"], s["response"], "x", color="black", ms=5, label="SF")
        ax.set_title(f"{k['name']}\nR² {k['r2']:.5f} | écart max {k['ecart_max_pct']:.1f} %", fontsize=8)
        ax.set_xlabel("ppm", fontsize=7); ax.set_ylabel("aire / aire ISTD", fontsize=7)
        ax.tick_params(labelsize=7); ax.legend(fontsize=6)
    for ax in list(axes.flat)[len(ok):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()

