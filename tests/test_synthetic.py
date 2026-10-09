"""Validation du traitement des pics sur des chromatogrammes synthétiques.

Pics EMG (gaussienne x exponentielle) pour reproduire la traînée à droite, bruit
blanc, dérive de baseline, m/z partagés et coélution partielle. La vérité terrain
(apex, aire, ratio) est connue par construction.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from scipy.stats import exponnorm

from reference import compute_reference

rng = np.random.default_rng(14)
DT = 0.004                       # pas d'échantillonnage (min)
T = np.arange(4.0, 19.0, DT)
SIG_W, TAU = 0.008, 0.05         # largeur gaussienne et constante de traînée (min)


def emg(loc, area):
    return area * exponnorm.pdf(T, TAU / SIG_W, loc=loc, scale=SIG_W)


# (nom, m/z quant, m/z qual, position, aire quant, ratio qual/quant vrai)
SPEC = [
    ("Isolé à traînée",     128, 127,  4.60, 900.0, 0.12),
    ("Paire séparée 1",     178, 176,  8.40, 700.0, 0.19),
    ("Paire séparée 2",     178, 176,  8.52, 650.0, 0.19),
    ("Coélution b",         252, 253, 13.30, 500.0, 0.22),
    ("Coélution k",         252, 253, 13.37, 480.0, 0.22),
    ("Troisième même m/z",  252, 253, 13.95, 450.0, 0.21),
    ("Faible signal",       276, 138, 16.20,  40.0, 0.25),
]


def build(noise=1.0):
    clean, traces = {}, {}
    for name, mq, ml, loc, area, ratio in SPEC:
        clean.setdefault(mq, np.zeros_like(T))
        clean.setdefault(ml, np.zeros_like(T))
        clean[mq] += emg(loc, area)
        clean[ml] += emg(loc, area * ratio)
    for mz, y in clean.items():
        drift = 5.0 + 0.4 * (T - T[0])                    # baseline lentement montante
        traces[mz] = (T, y + drift + rng.normal(0, noise, T.size))
    return traces


compounds = pd.DataFrame(
    [(n, mq, ml, i) for i, (n, mq, ml, *_ ) in enumerate(SPEC)],
    columns=["name", "mz_quant", "mz_qual", "elution_rank"])

ref = compute_reference(build(), compounds)

truth = {}
for name, mq, ml, loc, area, ratio in SPEC:
    y = emg(loc, area)
    truth[name] = (T[np.argmax(y)], area, ratio)

rows, ok = [], True
for _, r in ref.iterrows():
    rt_true, area_true, ratio_true = truth[r["name"]]
    e_rt = 100 * (r.rt_ref - rt_true) / rt_true
    e_area = 100 * (r.area_quant - area_true) / area_true
    e_ratio = 100 * (r.ratio_ref - ratio_true) / ratio_true
    rows.append((r["name"], r.rt_ref, e_rt, e_area, e_ratio, r.q_stop, r.warning))
    ok &= abs(e_rt) < 0.2 and abs(e_ratio) < 23
    if r["name"] in ("Isolé à traînée", "Troisième même m/z", "Faible signal"):
        ok &= abs(e_area) < 5          # pics isolés : l'aire doit être retrouvée
    elif abs(e_area) >= 5:
        # On ne prétend pas retrouver les aires coéluées : elles doivent être signalées.
        ok &= "séparation des aires à valider" in r.warning
print(pd.DataFrame(rows, columns=["composé", "rt_ref", "err_rt_%", "err_aire_%",
                                  "err_ratio_%", "arrêts", "alerte"])
      .to_string(index=False, float_format=lambda v: f"{v:.3f}"))

# Répétabilité : 200 tirages de bruit, dispersion du ratio et du temps de rétention.
rt, ra = [], []
for _ in range(200):
    r = compute_reference(build(), compounds).set_index("name")
    rt.append(r.rt_ref); ra.append(r.ratio_ref)
rt, ra = pd.concat(rt, axis=1), pd.concat(ra, axis=1)
rep = pd.DataFrame({"cv_rt_%": 100 * rt.std(axis=1) / rt.mean(axis=1),
                    "cv_ratio_%": 100 * ra.std(axis=1) / ra.mean(axis=1),
                    "biais_ratio_%": [100 * (ra.loc[n].mean() - truth[n][2]) / truth[n][2]
                                      for n in ra.index]})
print("\nRépétabilité sur 200 tirages de bruit :")
print(rep.to_string(float_format=lambda v: f"{v:.3f}"))
ok &= bool((rep["cv_rt_%"] < 0.2).all() and (rep["biais_ratio_%"].abs() < 23).all())
print("\nRESULTAT :", "OK (aires isolées et signalement des limites ; coélutions non validées)" if ok else "ECHEC")
sys.exit(0 if ok else 1)

