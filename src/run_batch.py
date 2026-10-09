"""Applique la référence GAM-6 aux 8 GAM et aux SF d'un batch.

Usage :
    python src/run_reference.py data/raw      # d'abord : la référence GAM-6
    python src/run_batch.py data/raw          # ensuite : tous les GAM et SF

Pour chaque échantillon et chaque composé : pic le plus proche du TR de
référence, contrôle TR (±0,2 %) et ratio qualifiant/quantifiant (±23 %), aires
quantifiant / qualifiant, aire de l'étalon interne et réponse normalisée.

Sorties dans outputs/<batch>/ :
    peaks.csv           une ligne par échantillon x composé
    conformite.csv      tableau composé x échantillon (OK / à vérifier / ...)
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from agilent_ms import export_sample                  # noqa: E402
from alerts import EDGE_MAX_RATIO, VALLEY_MAX_RATIO  # noqa: E402
from measure import measure_sample  # noqa: E402
from reference import default_method, load_compounds, load_sample  # noqa: E402

GAM_PPM = {1: 0.025, 2: 0.05, 3: 0.1, 4: 0.2, 5: 0.5, 6: 1.0, 7: 2.0, 8: 5.0}
SF_PPM = 1.0   # nominal de travail, à confirmer avec le cahier des charges


def classify(name: str):
    m = re.search(r"(GAM|SF)-.*-(\d+)$", name, re.I)
    if not m:
        return None
    kind, num = m.group(1).upper(), int(m.group(2))
    return kind, num, (GAM_PPM.get(num) if kind == "GAM" else SF_PPM)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", type=Path, help="dossier contenant le batch (dossiers .D)")
    ap.add_argument("--compounds", type=Path, default=default_method(ROOT),
                    help="metadata.xlsx de la méthode (défaut : data/metadata.xlsx ou ./metadata.xlsx)")
    ap.add_argument("--out", type=Path, default=ROOT / "outputs")
    ap.add_argument("--valley-max", type=float, default=None,
                    help="seuil vallée/hauteur (défaut : celui enregistré avec la référence)")
    ap.add_argument("--edge-max", type=float, default=None,
                    help="seuil de signal au bord SIM (défaut : celui enregistré avec la référence)")
    args = ap.parse_args()

    samples = sorted(p for p in args.source.resolve().rglob("*")
                     if p.is_dir() and p.suffix.lower() == ".d" and classify(p.stem))
    if not samples:
        sys.exit(f"Aucun dossier GAM ou SF .D trouvé sous {args.source}")
    parents = {p.parent for p in samples}
    if len(parents) != 1:
        sys.exit("Plusieurs batchs trouvés : indique un seul dossier batch pour éviter de mélanger les références")
    batch = samples[0].parent.name
    out = args.out / batch
    ref_path = out / "references_gam6.csv"
    if not ref_path.exists():
        sys.exit(f"{ref_path} introuvable : lance d'abord  python src/run_reference.py {args.source}")
    ref = pd.read_csv(ref_path)
    compounds = load_compounds(args.compounds)
    # Mêmes seuils que la référence, sauf choix explicite (et alors signalé).
    recorded = lambda col, default: float(ref[col].iloc[0]) if col in ref and ref[col].notna().any() else default
    valley_max = args.valley_max if args.valley_max is not None else recorded("valley_max", VALLEY_MAX_RATIO)
    edge_max = args.edge_max if args.edge_max is not None else recorded("edge_max", EDGE_MAX_RATIO)
    if (valley_max, edge_max) != (recorded("valley_max", valley_max), recorded("edge_max", edge_max)):
        print("ATTENTION : seuils différents de ceux de la référence GAM-6")
    print(f"Seuils expérimentaux : vallée <= {100 * valley_max:.0f} %, bord SIM <= {100 * edge_max:.1f} %")
    print(f"Méthode : {args.compounds.name} ({len(compounds)} composés)")

    tables = []
    for d in samples:
        kind, num, ppm = classify(d.stem)
        ion_dir = out / "ions" / d.stem
        if not any(ion_dir.glob("mz_*.csv")):
            export_sample(d, out / "ions")
        res = measure_sample(load_sample(ion_dir), ref, compounds, sample_type=kind,
                             valley_max=valley_max, edge_max=edge_max)
        res.insert(0, "sample", d.stem)
        res.insert(1, "type", kind)
        res.insert(2, "level", num)
        res.insert(3, "conc_nominal_ppm", ppm)
        tables.append(res)
        n_ok = (res["status"] == "OK").sum()
        print(f"{d.stem:28s} {kind:3s} {num}  -> {n_ok}/{len(res)} composés conformes")

    peaks = pd.concat(tables, ignore_index=True)
    peaks.to_csv(out / "peaks.csv", index=False, float_format="%.6g")
    conf = peaks.pivot_table(index="name", columns="sample", values="status", aggfunc="first")
    order = compounds["name"]
    conf = conf.reindex([n for n in order if n in conf.index])
    conf.to_csv(out / "conformite.csv")

    bad = peaks[peaks["status"] != "OK"]
    if len(bad):
        print("\nÀ vérifier :")
        cols = ["sample", "name", "status", "d_rt_pct", "d_ratio_pct", "blocking_alerts"]
        with pd.option_context("display.width", 220, "display.max_colwidth", 70):
            print(bad[[c for c in cols if c in bad]].to_string(index=False, float_format=lambda v: f"{v:+.2f}"))
    print(f"\nRésultats : {out / 'peaks.csv'}\nConformité : {out / 'conformite.csv'}")


if __name__ == "__main__":
    main()
