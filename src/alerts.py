"""Alertes de traitement : bloquantes ou informatives.

Une alerte **bloquante** empêche l'usage automatique d'une mesure (normalisation
par l'ISTD, entrée en calibration, verdict SF). Une alerte **informative** reste
visible dans les résultats mais ne bloque pas : elle décrit une situation
mesurée et jugée acceptable selon un critère explicite.

Par défaut, une situation est bloquante. Elle ne devient informative que si un
critère chiffré, configurable et documenté le justifie :

- vallée entre deux pics : informative si le rapport vallée/hauteur est
  inférieur ou égal à VALLEY_MAX_RATIO, bloquante sinon ;
- borne au bord d'une fenêtre SIM : informative si le signal résiduel au bord
  est inférieur ou égal à EDGE_MAX_RATIO de la hauteur, bloquante sinon ;
- borne arrêtée en « ligne de base » mais restée au-dessus de VALLEY_MAX_RATIO
  de la hauteur : bloquante (épaulement ou chevauchement non résolu).

Les deux seuils sont expérimentaux : ce sont des hypothèses de travail à tester
et à confirmer avec l'encadrement, pas une preuve que les aires sont justes.
Toutes les autres alertes (attribution ambiguë, pic réaffecté, ordre d'élution,
qualifiant absent ou décalé, double pic, filtrage des pics parasites, largeur
maximale) restent bloquantes.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

VALLEY_MAX_RATIO = 0.10   # expérimental, à confirmer : vallée <= 10 % de la hauteur
EDGE_MAX_RATIO = 0.01     # expérimental, à confirmer : signal au bord SIM <= 1 % de la hauteur


@dataclass(frozen=True)
class Alert:
    message: str
    blocking: bool = True


def peak_bound_alerts(peak, what: str = "", valley_max: float = VALLEY_MAX_RATIO,
                      edge_max: float = EDGE_MAX_RATIO) -> list[Alert]:
    """Alertes liées aux bornes d'un pic (vallée, bord SIM, largeur maximale)."""
    tag = f"{what} : " if what else ""
    out: list[Alert] = []
    if "limite" in (peak.left_stop, peak.right_stop):
        out.append(Alert(f"{tag}borne non trouvée (largeur maximale atteinte)"))
    v = peak.bound_ratio("vallee")
    if "vallee" in (peak.left_stop, peak.right_stop):
        if np.isfinite(v) and v <= valley_max:
            out.append(Alert(f"{tag}pic voisin séparé : vallée {100 * v:.1f} % de la hauteur "
                             f"(<= {100 * valley_max:.0f} %)", blocking=False))
        else:
            shown = "non mesurable" if not np.isfinite(v) else f"{100 * v:.1f} %"
            out.append(Alert(f"{tag}chevauchement : vallée {shown} de la hauteur "
                             f"(> {100 * valley_max:.0f} %), séparation des aires non validée"))
    # Une borne arrêtée en « ligne de base » doit être revenue près du fond. Si elle
    # reste au-dessus du seuil, le pic est sans doute un épaulement d'un voisin
    # (cas observé en simulation pour des pics très proches et un bruit élevé).
    b = peak.bound_ratio("baseline")
    if np.isfinite(b) and b > valley_max:
        out.append(Alert(f"{tag}borne sans retour au fond : niveau {100 * b:.1f} % de la hauteur "
                         f"(> {100 * valley_max:.0f} %), épaulement ou chevauchement possible"))
    e = peak.bound_ratio("fin_signal")
    if "fin_signal" in (peak.left_stop, peak.right_stop):
        if np.isfinite(e) and e <= edge_max:
            out.append(Alert(f"{tag}borne au bord de la fenêtre SIM : signal résiduel "
                             f"{100 * e:.2f} % de la hauteur (<= {100 * edge_max:.1f} %)", blocking=False))
        else:
            shown = "non mesurable" if not np.isfinite(e) else f"{100 * e:.2f} %"
            out.append(Alert(f"{tag}borne au bord de la fenêtre SIM : signal résiduel {shown} "
                             f"de la hauteur (> {100 * edge_max:.1f} %), aire tronquée possible"))
    return out


def render(alerts: list[Alert]) -> dict[str, str]:
    """Colonnes de sortie : toutes les alertes, puis séparées par gravité."""
    def join(items):
        seen, msgs = set(), []
        for a in items:
            if a.message not in seen:
                seen.add(a.message)
                msgs.append(a.message)
        return " ; ".join(msgs)
    return {"warning": join(alerts),
            "blocking_alerts": join([a for a in alerts if a.blocking]),
            "info_alerts": join([a for a in alerts if not a.blocking])}


def blocking_from_reference(row) -> list[Alert]:
    """Alertes bloquantes héritées de la référence GAM-6.

    Seules les alertes bloquantes de la référence sont propagées aux autres
    échantillons : les alertes informatives (vallée bien séparée, bord SIM à
    signal faible) sont réévaluées sur chaque échantillon. Si le fichier de
    référence ne distingue pas encore les gravités (ancienne version), toute
    alerte est traitée comme bloquante, par prudence.
    """
    if "blocking_alerts" in row.index:
        text = row.get("blocking_alerts")
    else:
        text = row.get("warning")
    if not isinstance(text, str) or not text.strip():
        return []
    return [Alert(f"référence : {m.strip()}") for m in text.split(" ; ") if m.strip()]
