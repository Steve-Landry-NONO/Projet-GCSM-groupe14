"""Application d'analyse GC-MS : référence, gamme, calibration, contrôle SF.

Lancement (depuis la racine du dépôt, venv actif) :
    streamlit run app.py

L'application ne recalcule rien : elle lit les sorties des scripts
(src/run_reference.py, src/run_batch.py, src/run_calibration.py) dans outputs/<batch>/.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

HERE = Path(__file__).resolve().parent
RT_TOL, RATIO_TOL = 0.2, 23.0          # % (cahier des charges, d'après les notes)
INK, MUTED, ACCENT, ALERT = "#1f2937", "#9ca3af", "#2563eb", "#dc2626"

st.set_page_config(page_title="Analyse GC-MS", page_icon="🧪", layout="wide")


# --------------------------------------------------------------------------- #
# Données
# --------------------------------------------------------------------------- #
@st.cache_data
def read_csv(path: str) -> pd.DataFrame | None:
    p = Path(path)
    return pd.read_csv(p) if p.exists() else None


@st.cache_data
def read_trace(batch_dir: str, sample: str, mz) -> pd.DataFrame | None:
    if pd.isna(mz):
        return None
    p = Path(batch_dir) / "ions" / sample / f"mz_{int(mz)}.csv"
    return pd.read_csv(p) if p.exists() else None


def fig_base(title="", x="", y="", height=360):
    f = go.Figure()
    f.update_layout(title=dict(text=title, font=dict(size=14)), height=height,
                    margin=dict(l=10, r=10, t=40 if title else 10, b=10),
                    xaxis_title=x, yaxis_title=y, template="plotly_white",
                    legend=dict(orientation="h", y=1.02, x=1, xanchor="right", yanchor="bottom"))
    return f


out_root = Path(st.sidebar.text_input("Dossier des résultats", str(HERE / "outputs")))
batches = sorted(d.name for d in out_root.iterdir() if d.is_dir()) if out_root.exists() else []
if not batches:
    st.error(f"Aucun résultat dans {out_root}. Lance d'abord src/run_reference.py, src/run_batch.py et src/run_calibration.py.")
    st.stop()
batch = st.sidebar.selectbox("Batch", batches, index=len(batches) - 1)
bdir = out_root / batch

ref = read_csv(str(bdir / "references_gam6.csv"))
peaks = read_csv(str(bdir / "peaks.csv"))
cal = read_csv(str(bdir / "calibrations.csv"))
pts = read_csv(str(bdir / "calibration_points.csv"))
sf = read_csv(str(bdir / "sf_results.csv"))

if ref is None:
    st.error("references_gam6.csv introuvable : lance run_reference.py.")
    st.stop()

names = ref.dropna(subset=["rt_ref"])["name"].tolist()
types = dict(zip(ref["name"], ref.get("type", pd.Series(index=ref.index, dtype=str))))
compound = st.sidebar.selectbox("Composé", names, index=names.index("Naphtalene") if "Naphtalene" in names else 0)
r0 = ref.set_index("name").loc[compound]
st.sidebar.caption(f"Type : **{types.get(compound, '?')}** · m/z {int(r0.mz_quant)} / {int(r0.mz_qual)}"
                   + (f" · ISTD : {r0.istd}" if isinstance(r0.get("istd"), str) else ""))

st.title("Analyse GC-MS — 16 HAP")
st.caption(f"Batch **{batch}** · référence **{r0['sample'] if 'sample' in r0 else 'GAM-6'}**")

tabs = st.tabs(["Vue d'ensemble", "Référence GAM-6", "Gamme GAM", "Calibration", "Contrôle SF", "BLPC"])

# --------------------------------------------------------------------------- #
# 0. Vue d'ensemble
# --------------------------------------------------------------------------- #
with tabs[0]:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Référence : composés trouvés", f"{ref['rt_ref'].notna().sum()} / {len(ref)}",
              help="Temps de rétention et ratio de référence trouvés sur GAM-6")
    if peaks is not None:
        g = peaks[peaks["type"] == "GAM"]
        c2.metric("Gamme GAM : mesures conformes", f"{(g['status'] == 'OK').sum()} / {len(g)}",
                  help="TR ±0,2 % et ratio ±23 % par rapport à la référence")
    else:
        c2.metric("Gamme GAM", "à lancer")
    c3.metric(f"Calibration : R² minimal ({len(cal)} HAP)" if cal is not None else "Calibration",
              f"{cal['r2'].min():.4f}" if cal is not None else "à lancer")
    if sf is not None and len(sf):
        c4.metric("Contrôle SF : résultats PASS", f"{sf['sf_status'].str.startswith('PASS').sum()} / {len(sf)}",
                  help="Tolérance fixée dans run_calibration.py (provisoire)")
    else:
        c4.metric("Contrôle SF", "à lancer")

    if peaks is not None:
        st.subheader("Conformité par composé et par échantillon")
        st.caption("Vert : TR et ratio dans les tolérances. Orange : à vérifier. Gris : non détecté.")
        piv = peaks.pivot_table(index="name", columns="sample", values="status", aggfunc="first")
        piv = piv.reindex([n for n in names if n in piv.index])
        piv.columns = [c.split("-", 1)[-1] for c in piv.columns]
        code = piv.apply(lambda col: col.map({"OK": 2, "à vérifier": 1})).fillna(0).astype(float)
        f = go.Figure(go.Heatmap(z=code.values, x=list(piv.columns), y=list(piv.index), zmin=0, zmax=2,
                                 colorscale=[[0, "#e5e7eb"], [0.5, "#f59e0b"], [1, "#16a34a"]],
                                 showscale=False, text=piv.values, hovertemplate="%{y}<br>%{x}<br>%{text}<extra></extra>",
                                 xgap=2, ygap=2))
        f.update_layout(height=40 + 22 * len(piv), margin=dict(l=10, r=10, t=10, b=10),
                        yaxis=dict(autorange="reversed"), template="plotly_white")
        st.plotly_chart(f, width="stretch")

# --------------------------------------------------------------------------- #
# 1. Référence
# --------------------------------------------------------------------------- #
with tabs[1]:
    st.markdown("Le temps de rétention de référence est l'apex du pic sur le m/z quantifiant ; "
                "les bornes viennent de la double descente depuis le sommet. "
                "Le ratio de référence = aire qualifiant / aire quantifiant.")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("TR de référence", f"{r0.rt_ref:.4f} min", help=f"Fenêtre ±0,2 % : {r0.rt_min:.4f} – {r0.rt_max:.4f}")
    m2.metric("Ratio de référence", f"{r0.ratio_ref:.4f}", help=f"Fenêtre ±23 % : {r0.ratio_min:.4f} – {r0.ratio_max:.4f}")
    m3.metric("Bornes", f"{r0.q_left:.3f} – {r0.q_right:.3f}")
    m4.metric("Signal / bruit", f"{r0.snr_quant:,.0f}".replace(",", " "))
    if isinstance(r0.get("warning"), str) and r0.warning:
        st.warning(r0.warning)

    sample = r0["sample"]
    tq, tl = read_trace(str(bdir), sample, r0.mz_quant), read_trace(str(bdir), sample, r0.mz_qual)
    if tq is None:
        st.info("Signaux par ion absents (dossier outputs/…/ions) : relance src/run_reference.py pour les regénérer.")
    else:
        width = max(0.25, r0.q_right - r0.q_left)
        lo, hi = st.slider("Fenêtre affichée (min)", float(tq.time_min.min()), float(tq.time_min.max()),
                           (float(r0.q_left - width / 2), float(r0.q_right + width / 2)), step=0.01)
        f = fig_base(x="temps (min)", y="intensité", height=420)
        f.add_vrect(x0=r0.q_left, x1=r0.q_right, fillcolor=MUTED, opacity=0.18, line_width=0,
                    annotation_text="bornes", annotation_position="top left")
        f.add_vrect(x0=r0.rt_min, x1=r0.rt_max, fillcolor=ACCENT, opacity=0.12, line_width=0)
        for tr, lab, dash, col in ((tq, f"quantifiant m/z {int(r0.mz_quant)}", "solid", INK),
                                   (tl, f"qualifiant m/z {int(r0.mz_qual)}", "dash", ACCENT)):
            if tr is not None:
                m = tr.time_min.between(lo, hi)
                f.add_scatter(x=tr.time_min[m], y=tr.intensity[m], name=lab, mode="lines",
                              line=dict(color=col, dash=dash, width=1.6))
        f.add_vline(x=r0.rt_ref, line=dict(color=INK, dash="dot", width=1))
        st.plotly_chart(f, width="stretch")
        st.caption("Gris : zone intégrée. Bleu clair : tolérance ±0,2 % autour du TR de référence.")

    st.subheader("Tous les composés")
    cols = ["elution_rank", "name", "type", "mz_quant", "mz_qual", "rt_ref", "rt_min", "rt_max",
            "q_left", "q_right", "ratio_ref", "ratio_min", "ratio_max", "warning"]
    st.dataframe(ref[[c for c in cols if c in ref]], hide_index=True, width="stretch")

# --------------------------------------------------------------------------- #
# 2. Gamme GAM
# --------------------------------------------------------------------------- #
with tabs[2]:
    if peaks is None:
        st.info("Lance src/run_batch.py pour mesurer les 8 GAM.")
    else:
        d = peaks[(peaks["name"] == compound)].sort_values(["type", "level"])
        g = d[d["type"] == "GAM"]
        left, right = st.columns(2)
        with left:
            f = fig_base("Écart de TR à la référence", "niveau GAM", "écart TR (%)")
            f.add_hrect(y0=-RT_TOL, y1=RT_TOL, fillcolor="#16a34a", opacity=0.08, line_width=0)
            f.add_scatter(x=g["level"], y=g["d_rt_pct"], mode="lines+markers", name="GAM",
                          line=dict(color=INK), marker=dict(color=np.where(g["rt_ok"] == True, INK, ALERT), size=9))
            s = d[d["type"] == "SF"]
            if len(s):
                f.add_scatter(x=[6] * len(s), y=s["d_rt_pct"], mode="markers", name="SF",
                              marker=dict(symbol="x", color=ACCENT, size=9))
            f.update_xaxes(tickvals=list(range(1, 9)))
            st.plotly_chart(f, width="stretch")
        with right:
            f = fig_base("Écart de ratio qualifiant/quantifiant", "niveau GAM", "écart ratio (%)")
            f.add_hrect(y0=-RATIO_TOL, y1=RATIO_TOL, fillcolor="#16a34a", opacity=0.08, line_width=0)
            f.add_scatter(x=g["level"], y=g["d_ratio_pct"], mode="lines+markers", name="GAM",
                          line=dict(color=INK), marker=dict(color=np.where(g["ratio_ok"] == True, INK, ALERT), size=9))
            f.update_xaxes(tickvals=list(range(1, 9)))
            st.plotly_chart(f, width="stretch")
        st.caption("Bande verte = tolérance. Point rouge = hors tolérance. "
                   "Niveaux 1 → 8 = 0,025 → 5 ppm ; la référence est le niveau 6 (1 ppm). "
                   "Les SF sont placées au niveau 6 (même concentration nominale).")

        st.subheader("Chromatogrammes superposés des 8 GAM")
        norm = st.toggle("Normaliser chaque courbe à son maximum", value=True)
        f = fig_base(x="temps (min)", y="intensité relative" if norm else "intensité", height=420)
        lo, hi = r0.q_left - 0.05, min(r0.q_right, r0.rt_ref + 0.25)
        shades = [f"rgba(37,99,235,{a:.2f})" for a in np.linspace(0.25, 1, 8)]
        found = False
        for (_, row), col in zip(g.iterrows(), shades):
            tr = read_trace(str(bdir), row["sample"], r0.mz_quant)
            if tr is None:
                continue
            found = True
            m = tr.time_min.between(lo, hi)
            y = tr.intensity[m] / (tr.intensity[m].max() if norm else 1)
            f.add_scatter(x=tr.time_min[m], y=y, mode="lines", line=dict(color=col, width=1.4),
                          name=f"GAM-{row['level']} ({row['conc_nominal_ppm']} ppm)")
        f.add_vline(x=r0.rt_ref, line=dict(color=INK, dash="dot", width=1))
        if found:
            st.plotly_chart(f, width="stretch")
        else:
            st.info("Signaux par ion des GAM absents (outputs/…/ions).")

        st.subheader("Stabilité des étalons internes")
        ist = peaks[peaks["type"].isin(["GAM", "SF"]) & peaks["name"].isin(ref.loc[ref.get("type") == "ISTD", "name"])]
        if len(ist):
            base = ist[ist["sample"].str.contains("GAM-.*-6$", regex=True)].set_index("name")["area_quant"]
            f = fig_base(y="aire / aire en GAM-6", height=320)
            for nm, dd in ist.groupby("name"):
                f.add_scatter(x=[s.split("-", 1)[-1] for s in dd["sample"]], y=dd["area_quant"] / base.get(nm, np.nan),
                              mode="lines+markers", name=nm)
            f.add_hline(y=1, line=dict(color=MUTED, dash="dot"))
            st.plotly_chart(f, width="stretch")
            st.caption("Une baisse commune à tous les étalons signale une perte de sensibilité de l'appareil ; "
                       "la normalisation aire / aire ISTD la corrige.")

        st.subheader("Mesures du composé")
        cols = ["sample", "level", "conc_nominal_ppm", "rt", "d_rt_pct", "ratio", "d_ratio_pct",
                "area_quant", "area_istd", "response", "status", "warning"]
        st.dataframe(d[[c for c in cols if c in d]], hide_index=True, width="stretch")

# --------------------------------------------------------------------------- #
# 3. Calibration
# --------------------------------------------------------------------------- #
with tabs[3]:
    if cal is None:
        st.info("Lance src/run_calibration.py.")
    elif compound not in set(cal["name"]):
        st.info(f"{compound} n'est pas calibré (étalon interne ou surrogat). Choisis un HAP cible.")
    else:
        k = cal.set_index("name").loc[compound]
        p = pts[pts["name"] == compound]
        c1, c2, c3 = st.columns(3)
        c1.metric("R²", f"{k.r2:.5f}")
        c2.metric("Écart max des points recalculés", f"{k.ecart_max_pct:.1f} %")
        c3.metric("Points de gamme", f"{int(k.n_points)}")
        st.markdown(f"**Modèle** (pondération {k.ponderation}) : réponse = "
                    f"`{k.a:.5f}·C² {k.b:+.5f}·C {k.c:+.5f}`, avec réponse = aire / aire {k.istd}")
        left, right = st.columns([3, 2])
        with left:
            log = st.toggle("Échelle logarithmique", value=False,
                            help="Utile pour voir les bas niveaux, tassés près de zéro en échelle linéaire")
            xs = np.geomspace(max(k.c_min_ppm / 2, 1e-3), k.c_max_ppm, 200) if log else np.linspace(0, k.c_max_ppm, 200)
            f = fig_base("Courbe de calibration", "concentration (ppm)", "aire / aire ISTD", height=420)
            f.add_scatter(x=xs, y=np.polyval([k.a, k.b, k.c], xs), mode="lines", name="ajustement",
                          line=dict(color=INK, width=1.5))
            f.add_scatter(x=p["conc_nominal_ppm"], y=p["response"], mode="markers", name="GAM",
                          marker=dict(color="white", line=dict(color=INK, width=1.5), size=9))
            if sf is not None and len(sf):
                s = sf[sf["name"] == compound]
                f.add_scatter(x=s["conc_final_ppm"], y=s["response"], mode="markers", name="SF (recalculées)",
                              marker=dict(symbol="x", color=ACCENT, size=10))
            if log:
                f.update_xaxes(type="log"); f.update_yaxes(type="log")
            st.plotly_chart(f, width="stretch")
        with right:
            f = fig_base("Écart des points recalculés", "niveau GAM", "écart (%)", height=420)
            f.add_bar(x=p["level"], y=p["ecart_pct"], marker_color=np.where(p["ecart_pct"].abs() > 15, ALERT, INK))
            f.update_xaxes(tickvals=list(range(1, 9)))
            st.plotly_chart(f, width="stretch")
        st.caption("Écart = (concentration recalculée par la courbe − concentration connue) / concentration connue. "
                   "Un R² proche de 1 ne suffit pas : les bas niveaux pèsent peu dans le R² et se lisent ici.")
        st.subheader("Tous les HAP")
        st.dataframe(cal, hide_index=True, width="stretch")

# --------------------------------------------------------------------------- #
# 4. Contrôle SF
# --------------------------------------------------------------------------- #
with tabs[4]:
    if sf is None or not len(sf):
        st.info("Lance src/run_calibration.py.")
    else:
        tol = st.slider("Tolérance d'acceptation (± %)", 5, 40, 20, 1,
                        help="Valeur provisoire : la tolérance officielle reste à confirmer avec le cahier des charges.")
        s = sf.copy()
        s["verdict"] = np.where(s["ecart_pct"].abs() <= tol, "PASS", "FAIL")
        n_pass = (s["verdict"] == "PASS").sum()
        c1, c2, c3 = st.columns(3)
        c1.metric("Résultats conformes", f"{n_pass}/{len(s)}")
        c2.metric("Écart médian", f"{s['ecart_pct'].median():+.1f} %")
        c3.metric("Nominal", f"{s['conc_nominal_ppm'].iloc[0]} ppm", help="Valeur de travail, à confirmer")
        piv = s.pivot_table(index="name", columns="sample", values="ecart_pct")
        piv = piv.reindex([n for n in names if n in piv.index])
        piv.columns = [c.split("-", 1)[-1] for c in piv.columns]
        lim = max(tol * 1.5, float(np.nanmax(np.abs(piv.values))))
        f = go.Figure(go.Heatmap(z=piv.values, x=list(piv.columns), y=list(piv.index), zmid=0, zmin=-lim, zmax=lim,
                                 colorscale="RdBu",
                                 text=[[("—" if np.isnan(v) else f"{v:+.1f} %") for v in row] for row in piv.values],
                                 texttemplate="%{text}",
                                 hovertemplate="%{y}<br>%{x}<br>écart %{z:.1f} %<extra></extra>",
                                 colorbar=dict(title="écart %"), xgap=2, ygap=2))
        f.update_layout(height=60 + 24 * len(piv), margin=dict(l=10, r=10, t=10, b=10),
                        yaxis=dict(autorange="reversed"), template="plotly_white")
        st.plotly_chart(f, width="stretch")
        st.caption("Écart = (concentration retrouvée − nominal) / nominal. Un écart du même signe sur tous les HAP "
                   "suggère une cause commune (nominal, préparation, étalon interne) plutôt qu'un problème de pic.")
        st.dataframe(s, hide_index=True, width="stretch")

# --------------------------------------------------------------------------- #
# 5. BLPC
# --------------------------------------------------------------------------- #
with tabs[5]:
    st.info("Les BLPC (échantillons réels) sont hors du périmètre actuel du cahier des charges. "
            "La chaîne est prête à les traiter de la même manière : mesure contre la référence GAM-6, "
            "inversion de la calibration, puis application du facteur de dilution (multiplier) une seule fois.")
