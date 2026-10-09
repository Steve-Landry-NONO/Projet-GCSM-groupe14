"""Caractérise le critère vallée/hauteur sur une grille de cas synthétiques.

    python tests/caracteriser_vallee.py [--seeds 3] [--valley-max 0.10]

Affiche les tableaux repris dans docs/VALLEE_HAUTEUR.md. Les chiffres
dépendent des formes de pics simulées : ils décrivent le comportement de
l'algorithme, pas la justesse des aires sur les échantillons réels.
"""
import argparse
import itertools
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from alerts import peak_bound_alerts  # noqa: E402
from valley_cases import GRID, run_case  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--valley-max", type=float, default=0.10)
    a = ap.parse_args()

    rows = []
    for hr, sep, noise, base in itertools.product(*GRID.values()):
        for seed in range(a.seeds):
            rows += run_case(hr, sep, noise, base, seed, with_peaks=True)
    d = pd.DataFrame(rows)
    print(f"Cas simulés : {len(d)} pics ({len(d) // 2} couples)")
    print(f"Couples non séparés par la détection (comptés comme bloquants en amont) : "
          f"{(~d.detected).sum() // 2}")
    d = d[d.detected].copy()
    d["blocked"] = [any(x.blocking for x in peak_bound_alerts(p, valley_max=a.valley_max))
                    for p in d.pop("peak_obj")]
    d["err_neighbour"] = d.err_area - d.err_iso      # part de l'erreur due au voisin

    print("\nPrécision du rapport mesuré (écart au cas sans bruit, en points de %) :")
    t = (100 * (d.v_est - d.v_ref).abs()).groupby([d.noise, d.base]).agg(["median", "max"])
    print(t.round(1).to_string())

    print("\nDécision et erreur d'aire, par niveau de bruit :")
    out = []
    for nz, s in d.groupby("noise"):
        free = s[~s.blocked]
        out.append({
            "bruit (% hauteur max)": 100 * nz,
            "pics mesurés": len(s),
            "non bloqués": len(free),
            "creux réel > 15 % non bloqués": int((free.v_ref > 0.15).sum()),
            "|err aire| max non bloqués (%)": round(free.err_area.abs().max(), 1),
            "|err due au voisin| médiane (%)": round(free.err_neighbour.abs().median(), 1),
            "|err due au voisin| max (%)": round(free.err_neighbour.abs().max(), 1),
            "non bloqués avec |err voisin| > 5 %": int((free.err_neighbour.abs() > 5).sum()),
            "|err| témoin isolé médiane (%)": round(s.err_iso.abs().median(), 1),
        })
    print(pd.DataFrame(out).to_string(index=False))

    print("\nBruit 0,1 %, ligne plate, tirage 0 — rapport et erreur d'aire :")
    s = d[(d.noise == 0.001) & (d.base == "plate") & (d.seed == 0)]
    print(s[["height_ratio", "sep", "peak", "v_ref", "v_est", "err_area", "blocked"]]
          .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
