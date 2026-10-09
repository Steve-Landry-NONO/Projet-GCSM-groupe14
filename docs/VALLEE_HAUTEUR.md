# Rapport vallée/hauteur et classement des alertes

*9 octobre 2026 — critère expérimental, à confirmer avec l'encadrement.*

## Pourquoi ce critère

Quand deux composés sortent l'un près de l'autre sur le même ion, la descente depuis l'apex s'arrête au creux entre les deux pics (une « vallée ») : l'aire est séparée par une coupure verticale. Cette coupure attribue à chaque pic une partie de la traînée de l'autre. L'erreur est faible si le creux descend presque jusqu'au fond, forte s'il reste haut.

Avant ce critère, toute vallée bloquait la mesure. Sur le rejeu du batch 20251103, aucune réponse stricte n'était alors calculable, alors que les creux observés vont de moins de 1 % à environ 50 % de la hauteur selon les paires. Le rapport vallée/hauteur mesure cette profondeur pour distinguer les cas.

## Définition du calcul

Pour une borne b d'un pic :

```text
v = (L_b - F) / (H - F)        borné à [0, 1]
```

| Terme | Calcul (`src/peaks.py`) |
| --- | --- |
| `L_b` | Niveau du signal brut à la borne : médiane des 5 points centrés sur la borne. |
| `H` | Niveau du signal lissé (Savitzky-Golay, 7 points) à l'apex. |
| `F` | Fond : le plus bas entre les niveaux des deux bornes et le 10e percentile du signal brut à ±0,5 min autour de l'apex. |

- `v = 0` : le signal revient au fond entre les deux pics.
- `v = 1` : aucun creux, le pic est un épaulement.
- Le rapport est calculé **par rapport à la hauteur du pic considéré**. Pour un petit pic accolé à un grand, il est donc plus élevé pour le petit, qui est aussi celui dont l'aire est la plus affectée.
- Le 10e percentile évite de prendre un niveau de vallée comme fond quand un pic est encadré par deux voisins. Prendre le minimum rend le rapport plutôt surestimé, donc prudent.
- Le même calcul donne le **signal résiduel au bord d'une fenêtre SIM**, quand l'acquisition de l'ion s'arrête avant le retour au fond.

Les valeurs sont écrites dans `references_gam6.csv` et `peaks.csv` (colonnes `valley_ratio`, `edge_ratio`), avec les seuils utilisés (`valley_max`, `edge_max`).

## Alertes bloquantes et informatives

Une alerte bloquante empêche la normalisation par l'étalon interne, l'entrée en calibration et tout verdict SF. Une alerte informative reste affichée mais ne bloque pas. **Par défaut, toute alerte est bloquante.** Elle ne devient informative que sur un critère chiffré et configurable.

| Situation | Informative si | Sinon |
| --- | --- | --- |
| Borne à une vallée | `v <= valley_max` (défaut 0,10) | Bloquante (chevauchement, ou rapport non mesurable) |
| Borne au bord d'une fenêtre SIM | signal résiduel `<= edge_max` (défaut 0,01) | Bloquante (aire tronquée possible) |
| Borne « ligne de base » restée haute | jamais | Bloquante si `v > valley_max` (épaulement) |
| Attribution ambiguë (« autre pic proche en intensité ») | jamais | Bloquante |
| Pic réaffecté par l'ordre d'élution, ordre non respecté | jamais | Bloquante |
| Qualifiant absent, décalé, faible ou fenêtre absente | jamais | Bloquante |
| `double_peak`, `median_despike` | jamais | Bloquante (règle non implémentée) |
| Largeur maximale atteinte, apex attribué deux fois, ISTD non validé | jamais | Bloquante |

Depuis la référence GAM-6, **seules les alertes bloquantes sont propagées** aux autres échantillons. Les alertes informatives sont recalculées sur chaque échantillon, puisque le creux dépend des concentrations. Une référence produite par une version antérieure, sans gravité, est traitée comme entièrement bloquante.

Seuils en ligne de commande :

```bash
python src/run_reference.py data/raw/20251103 --valley-max 0.10 --edge-max 0.01
python src/run_batch.py data/raw/20251103      # reprend les seuils enregistrés avec la référence
```

## Caractérisation sur cas synthétiques

`tests/caracteriser_vallee.py` simule des couples de pics gaussiens à traînée exponentielle, au pas réel du batch (0,0037 min). Il fait varier :
- les proportions de hauteur : 1, 0,2, 5 et 0,05 ;
- les écarts entre apex : 0,045 à 0,18 min ;
- le bruit : 0,1 %, 1 % et 3 % de la plus grande hauteur ;
- la ligne de base : plate, en pente, décalée.

Chaque combinaison est tirée 3 fois. Un témoin « pic seul » sépare l'erreur due au voisin de l'erreur due au bruit.

Sur ces 756 couples, 333 ne sont pas séparés par la détection. Ce sont surtout les pics très proches, les très petits pics et le bruit élevé. Ils sont bloqués en amont par l'alerte « n pics pour n composés ».

**Précision du rapport mesuré** (écart au même cas sans bruit, en points de pourcentage) :

| Bruit | Ligne plate (médiane / max) | En pente | Décalée |
| --- | --- | --- | --- |
| 0,1 % | 0,1 / 2,1 | 0,9 / 10,6 | 0,4 / 6,1 |
| 1 % | 1,2 / 8,6 | 1,5 / 10,7 | 1,3 / 9,1 |
| 3 % | 2,6 / 5,0 | 2,8 / 5,2 | 2,6 / 5,0 |

**Décision et erreur d'aire des pics non bloqués :**

| Bruit | Pics mesurés | Non bloqués | Creux réel > 15 % non bloqués | Erreur d'aire max | Part due au voisin (médiane / max) | Non bloqués avec erreur due au voisin > 5 % |
| --- | --- | --- | --- | --- | --- | --- |
| 0,1 % | 414 | 183 | 0 | 6,3 % | 1,6 % / 8,1 % | 7 |
| 1 % | 306 | 140 | 0 | 9,1 % | 2,6 % / 8,1 % | 26 |
| 3 % | 126 | 46 | 0 | 16,1 % | 3,7 % / 13,7 % | 15 |

## Ce que ces résultats montrent, et ce qu'ils ne montrent pas

- **Aucun chevauchement marqué ne passe.** Aucun creux réel supérieur à 15 % de la hauteur n'est classé informatif, quels que soient le bruit et la ligne de base (test `test_no_marked_overlap_becomes_informative`). Les épaulements sans creux mesurable sont bloqués par la règle « borne sans retour au fond ».
- **Un pic non bloqué n'a pas une aire garantie.** Sous le seuil de 10 %, la coupure à la vallée laisse une erreur due au voisin d'environ 2 % en médiane, mais jusqu'à 8 % à faible bruit et 14 % à bruit élevé dans ces simulations. **Le seuil de 10 % est une hypothèse de travail, pas une preuve de justesse des aires.**
- **La précision du rapport baisse avec la pente et le bruit.** Près du seuil, des écarts de 5 à 10 points peuvent faire basculer la décision d'un échantillon à l'autre.
- **Le bruit seul dégrade aussi l'aire.** Le témoin isolé perd en médiane 0,7 % (bruit 0,1 %) à 6,7 % (bruit 3 %) de son aire, parce que la descente s'arrête plus tôt dans la traînée bruitée. Cette limite concerne tous les pics, voisins ou non.
- **Les simulations ne remplacent pas les données.** Elles utilisent une seule forme de pic (gaussienne à traînée). Les pics réels peuvent être plus asymétriques. Sur GAM-6, les creux rapportés sont d'environ 0,8 à 1 % pour Phénanthrène/Anthracène, 1,5 à 1,9 % pour Benzo(a)anthracène/Chrysène et 47 à 52 % pour les Benzo(b)/(k)fluoranthènes. C'est le chevauchement le plus marqué parmi les paires examinées, sans que cela prouve qu'il n'y en ait pas d'autre dans d'autres échantillons ou à d'autres concentrations.

## À confirmer avec l'encadrement

1. Le principe d'un seuil vallée/hauteur, et sa valeur.
2. La méthode de séparation attendue pour une coélution marquée : coupure verticale, tangente, ajustement de profils, ou bornes MassHunter.
3. Le seuil de signal résiduel acceptable au bord d'une fenêtre SIM.
4. La comparaison bornes à bornes avec les QuantReports MassHunter, qui reste à faire avant toute acceptation quantitative.
