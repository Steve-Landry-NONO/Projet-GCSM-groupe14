"""Lecture des fichiers Agilent `data.ms` (dossiers .D) et export des ions en CSV.

Format (ChemStation MSD, big-endian) : pour chaque scan, un temps en ms/60 et une
liste de paires (m/z x 20, abondance codée mantisse 14 bits x 8^exposant 2 bits).
En mode SIM, un scan ne contient que les ions du groupe actif : un ion absent
d'un scan n'a pas été mesuré, il ne vaut pas zéro. L'export ne garde donc, pour
chaque ion, que les instants où il a réellement été acquis.

Usage :
    python src/agilent_ms.py chemin/vers/ECHANTILLON.D  dossier_sortie  [m/z ...]
    -> dossier_sortie/ECHANTILLON/mz_128.csv  (colonnes time_min, intensity)
"""
from __future__ import annotations

import struct
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def read_data_ms(path: str | Path) -> pd.DataFrame:
    """Retourne une table longue : scan, time_min, mz (mesuré), intensity."""
    raw = Path(path).read_bytes()
    if struct.unpack_from(">I", raw, 0)[0] != 0x01320000:
        raise ValueError(f"{path} : en-tête data.ms Agilent non reconnu")
    kind = raw[5:5 + raw[4]].decode("latin-1").strip()
    if kind == "MSD Spectral File":                      # variante LC
        n_scans = struct.unpack_from(">I", raw, 0x116)[0]
    else:                                                # GC/MS
        n_scans = struct.unpack_from("<H", raw, 0x142)[0]
    pos = struct.unpack_from(">H", raw, 0x10A)[0] * 2 - 2

    scans, times, mzs, ints = [], [], [], []
    for i in range(n_scans):
        t_raw = struct.unpack_from(">I", raw, pos + 2)[0]
        n = struct.unpack_from(">H", raw, pos + 12)[0]
        pairs = np.frombuffer(raw, dtype=">u2", count=2 * n, offset=pos + 18).reshape(n, 2)
        enc = pairs[:, 1].astype(np.int64)
        scans.append(np.full(n, i))
        times.append(np.full(n, t_raw / 60000.0))
        mzs.append(pairs[:, 0] / 20.0)
        ints.append((enc & 0x3FFF) * 8 ** (enc >> 14))
        pos += 18 + 4 * n + 10
    return pd.DataFrame({"scan": np.concatenate(scans), "time_min": np.concatenate(times),
                         "mz": np.concatenate(mzs), "intensity": np.concatenate(ints)})


def ion_traces(table: pd.DataFrame) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """{m/z nominal: (time_min, intensity)}, limité aux scans où l'ion est acquis."""
    table = table.assign(mz_nom=np.rint(table["mz"]).astype(int))
    g = table.groupby(["mz_nom", "scan"], sort=True).agg(
        time_min=("time_min", "first"), intensity=("intensity", "sum")).reset_index()
    return {int(mz): (d["time_min"].to_numpy(), d["intensity"].to_numpy(float))
            for mz, d in g.groupby("mz_nom")}


def export_sample(d_folder: str | Path, out_dir: str | Path, mz_list=None) -> Path:
    """Écrit un mz_XXX.csv par ion pour un dossier .D ; retourne le dossier créé."""
    d_folder = Path(d_folder)
    ms = next((p for p in d_folder.iterdir() if p.name.lower() == "data.ms"), None)
    if ms is None:
        raise FileNotFoundError(f"data.ms introuvable dans {d_folder}")
    traces = ion_traces(read_data_ms(ms))
    dest = Path(out_dir) / d_folder.name.removesuffix(".D").removesuffix(".d")
    dest.mkdir(parents=True, exist_ok=True)
    for mz, (t, y) in traces.items():
        if mz_list and mz not in mz_list:
            continue
        pd.DataFrame({"time_min": t, "intensity": y}).to_csv(
            dest / f"mz_{mz}.csv", index=False, float_format="%.6f")
    return dest


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    wanted = {int(m) for m in sys.argv[3:]} or None
    print(export_sample(sys.argv[1], sys.argv[2], wanted))

