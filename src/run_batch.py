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
    args = ap.parse_args()

    samples = sorted(p for p in args.source.resolve().rglob("*")
                     if p.is_dir() and p.suffix.lower() == ".d" and classify(p.stem))
    if not samples:
        sys.exit(f"Aucun dossier GAM ou SF .D trouvé sous {args.source}")
    batch = samples[0].parent.name
    out = args.out / batch
    ref_path = out / "references_gam6.csv"
    if not ref_path.exists():
        sys.exit(f"{ref_path} introuvable : lance d'abord  python src/run_reference.py {args.source}")
    ref = pd.read_csv(ref_path)
    compounds = load_compounds(args.compounds)
    print(f"Méthode : {args.compounds.name} ({len(compounds)} composés)")

    tables = []
    for d in samples:
        kind, num, ppm = classify(d.stem)
        ion_dir = out / "ions" / d.stem
        if not any(ion_dir.glob("mz_*.csv")):
            export_sample(d, out / "ions")
        res = measure_sample(load_sample(ion_dir), ref, compounds)
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
    order = ref.dropna(subset=["rt_ref"])["name"]
    conf = conf.reindex([n for n in order if n in conf.index])
    conf.to_csv(out / "conformite.csv")

    bad = peaks[peaks["status"] != "OK"]
    if len(bad):
        print("\nÀ vérifier :")
        cols = ["sample", "name", "status", "d_rt_pct", "d_ratio_pct", "warning"]
        with pd.option_context("display.width", 220, "display.max_colwidth", 70):
            print(bad[[c for c in cols if c in bad]].to_string(index=False, float_format=lambda v: f"{v:+.2f}"))
    print(f"\nRésultats : {out / 'peaks.csv'}\nConformité : {out / 'conformite.csv'}")


if __name__ == "__main__":
    main()
