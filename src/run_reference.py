"""Calcule la référence GAM-6 : temps de rétention et ratio qualifiant/quantifiant.

Usage (depuis le dossier du projet) :
    python src/run_reference.py data/raw/20251103
    python src/run_reference.py data/raw/20251103/20251103-GAM-25-385-6.D

Entrée : un dossier de batch (le GAM-6 y est trouvé automatiquement), un dossier
.D, ou un dossier de mz_XXX.csv déjà préparés.
Sorties dans outputs/<batch>/ :
    ions/<échantillon>/mz_XXX.csv   signaux par ion extraits du data.ms
    references_gam6.csv             une ligne par composé
    controle_gam6.png               pic quantifiant (plein) + qualifiant, bornes
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
from agilent_ms import export_sample                     # noqa: E402
from reference import compute_reference, default_method, load_compounds, load_sample  # noqa: E402


def find_gam6(folder: Path) -> Path:
    if folder.suffix.lower() == ".d" or any(folder.glob("mz_*.csv")):
        return folder
    hits = [p for p in folder.rglob("*") if p.is_dir()
            and re.search(r"GAM.*-6(\.D)?$", p.name, re.I)
            and (p.suffix.lower() == ".d" or any(p.glob("mz_*.csv")))]
    if len(hits) != 1:
        found = "\n  ".join(str(h) for h in hits) or "aucun"
        sys.exit(f"Impossible de choisir le GAM-6 dans {folder} ; candidats :\n  {found}\n"
                 "Donne directement le chemin du dossier .D.")
    return hits[0]


def plot(ref: pd.DataFrame, traces: dict, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ok = ref.dropna(subset=["rt_ref"])
    cols = 4
    if ok.empty:
        path.unlink(missing_ok=True)
        return
    rows = -(-len(ok) // cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 2.9 * rows), squeeze=False)
    for ax, (_, r) in zip(axes.flat, ok.iterrows()):
        pad = max(0.15, (r.q_right - r.q_left))
        lo, hi = r.q_left - pad / 2, r.q_right + pad
        for mz, style, lab in ((int(r.mz_quant), "-", "quant"), (r.mz_qual, "--", "qual")):
            if pd.isna(mz) or int(mz) not in traces:
                continue
            t, y = traces[int(mz)]
            m = (t >= lo) & (t <= hi)
            ax.plot(t[m], y[m], style, color="black", lw=0.9, label=f"{lab} m/z {int(mz)}")
        ax.axvspan(r.q_left, r.q_right, color="0.85", zorder=0)
        ax.axvline(r.rt_ref, color="black", lw=0.6, ls=":")
        title = f"{r['name']}\nTR {r.rt_ref:.4f} min | ratio {r.get('ratio_ref', float('nan')):.3f}"
        if isinstance(r.warning, str) and r.warning:
            title += "  (!)"
        ax.set_title(title, fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, loc="upper right")
    for ax in list(axes.flat)[len(ok):]:
        ax.axis("off")
    fig.suptitle("GAM-6 : zone grise = bornes d'intégration, pointillé = temps de rétention", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", type=Path, help="dossier batch, dossier .D ou dossier de mz_XXX.csv")
    ap.add_argument("--compounds", type=Path, default=default_method(ROOT),
                    help="metadata.xlsx de la méthode (défaut : data/metadata.xlsx ou ./metadata.xlsx)")
    ap.add_argument("--out", type=Path, default=ROOT / "outputs")
    ap.add_argument("--qual-bounds", choices=["quantifiant", "independant"], default="quantifiant")
    args = ap.parse_args()

    src = find_gam6(args.source.resolve())
    batch = src.parent.name if src.suffix.lower() == ".d" else args.source.resolve().name
    out = args.out / batch
    out.mkdir(parents=True, exist_ok=True)

    if src.suffix.lower() == ".d":
        print(f"Extraction des ions de {src.name} ...")
        ion_dir = export_sample(src, out / "ions")
    else:
        ion_dir = src
    traces = load_sample(ion_dir)
    print(f"{len(traces)} ions lus : {sorted(traces)}")

    compounds = load_compounds(args.compounds)
    print(f"Méthode : {args.compounds.name} ({len(compounds)} composés)")
    missing = sorted((set(compounds.mz_quant) | set(compounds.mz_qual.dropna().astype(int))) - set(traces))
    if missing:
        print(f"ATTENTION : m/z de la méthode absents des données : {missing}")

    ref = compute_reference(traces, compounds, qual_bounds=args.qual_bounds)
    ref.insert(0, "sample", ion_dir.name)
    extra = [c for c in ("type", "use", "istd") if c in compounds]
    ref = ref.merge(compounds[["name"] + extra], on="name", how="left")
    ref.to_csv(out / "references_gam6.csv", index=False, float_format="%.6g")
    plot(ref, traces, out / "controle_gam6.png")

    show = ["elution_rank", "name", "mz_quant", "mz_qual", "rt_ref", "q_left", "q_right",
            "ratio_ref", "warning"]
    with pd.option_context("display.width", 200, "display.max_colwidth", 60):
        print(ref[[c for c in show if c in ref]].to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"\nRésultats : {out / 'references_gam6.csv'}\nContrôle visuel : {out / 'controle_gam6.png'}")


if __name__ == "__main__":
    main()

