# Contexte du projet GC-MS — groupe 14

Ce fichier sert à reprendre le projet sans relire les conversations : il rassemble le contexte métier, l'état du travail, les décisions prises et les questions ouvertes. Il est tenu à jour à chaque étape importante, par les membres du groupe comme par les assistants (Claude, ChatGPT).

*Dernière mise à jour : 9 octobre 2026.*

## 1. Le sujet en bref

Un laboratoire mesure la concentration de 16 HAP (hydrocarbures aromatiques polycycliques, cancérigènes) dans des échantillons de terre ou d'enrobé, avec une GC-MS (chromatographie en phase gazeuse couplée à la spectrométrie de masse). Aujourd'hui, le chimiste repère les pics, fixe les bornes d'intégration et calcule les aires à la main dans MassHunter. **Le projet consiste à automatiser ce travail** de façon reproductible et explicable.

### La chaîne physique

1. L'échantillon de terre est dissous dans un solvant, l'hexane. On n'injecte qu'une petite partie de la solution (par exemple 1 ml dans une fiole de 10 ml, d'où un facteur de dilution).
2. **Chromatographie** : la solution traverse une colonne capillaire (silice, environ 20 m) placée dans un four. Les HAP adhèrent plus ou moins à la paroi et sortent à des temps différents : c'est le **temps de rétention**.
3. **Spectrométrie de masse** : chaque molécule qui sort est fragmentée et ionisée. Le détecteur compte les fragments selon leur rapport masse/charge (**m/z**). Chaque molécule a une empreinte de fragmentation propre (sa « loi de probabilité » de fragmentation).

### Deux invariants à retenir

- **Le temps de rétention absolu varie** avec l'état de la colonne (dépôts, usure) et, on l'a vérifié, avec la concentration. **L'ordre de sortie, lui, ne change pas.**
- Il faut séparer deux types de données : les **invariants** de la littérature (ions caractéristiques, ordre d'élution) et les **données d'expérience** (signaux mesurés, temps observés).

### Vocabulaire

| Terme | Sens |
| --- | --- |
| Ion quantifiant | m/z utilisé pour mesurer l'aire du composé |
| Ion qualifiant | second m/z, qui sert à confirmer l'identité du pic |
| Ratio qualifiant/quantifiant | aire qualifiant / aire quantifiant. Il est caractéristique d'une molécule et sert de contrôle d'identité. |
| ISTD (étalon interne) | molécule deutérée ajoutée en quantité fixe ; on normalise par son aire |
| Deutéré (ex. Naphtalène-D8) | molécule dont les 8 hydrogènes sont remplacés par du deutérium (hydrogène avec un neutron de plus). La chimie est quasi identique, mais la masse est plus élevée (m/z 136 au lieu de 128). |
| Surrogat | traceur ajouté pour suivre l'injection ou l'extraction (Pyrène-D10, Benzo(a)anthracène-D12) |
| Coélution | deux composés sortent presque en même temps, avec des pics qui se chevauchent |
| SIM | mode d'acquisition où seuls certains ions sont mesurés, par fenêtres de temps |

## 2. La méthode du chimiste, et donc le pipeline

| Étape | Échantillons | Rôle | Tolérance |
| --- | --- | --- | --- |
| **Référence** | GAM-6 (1 ppm) | Temps de rétention et ratio qualifiant/quantifiant de référence de chaque composé | TR ±0,2 % ; ratio ±23 % |
| **Gamme (GAM)** | 8 niveaux : 0,025 ; 0,05 ; 0,1 ; 0,2 ; 0,5 ; 1 ; 2 ; 5 ppm | Calibration quadratique par HAP : réponse = a·C² + b·C + c | — |
| **Contrôle (SF)** | solutions synthétiques à 1 ppm | Vérifier que la calibration retrouve la bonne concentration | à confirmer |
| **BLPC** | échantillons réels | Application finale | hors périmètre actuel |

La réponse est le rapport **aire quantifiant du composé / aire quantifiant de son ISTD**.

Pour chaque échantillon, on retient le pic le plus proche du temps de référence. Il est accepté si son TR est dans ±0,2 % et son ratio dans ±23 % de la référence.

Analogie donnée en cours : la gamme est le jeu d'entraînement, les SF le jeu de test, les BLPC l'application.

## 3. Les données

- **Batch 20251103** : 46 acquisitions (8 GAM, 4 SF, 6 hexane, 28 BLPC) au format Agilent `.D`.
  - Elles sont hors du repo (volumineuses, provenance instrumentale). Il faut les placer dans `data/raw/`.
  - Les signaux sont dans `data.ms` : 33 ions en SIM, avec un pas de temps d'environ 0,0037 min.
- **`metadata.xlsx`** : la méthode, et la source de vérité pour les composés. Il est hors du repo et va dans `data/metadata.xlsx`.
  - Feuille `compound_informations` : 23 composés (16 cibles, 4 ISTD, 3 surrogats), avec ions quantifiant et qualifiant, `double_peak`, `coelution_order` et `median_despike`.
  - **L'ordre des lignes est l'ordre d'élution**, ce qui est confirmé par les données.
  - Feuille `compound_deuterated` : l'ISTD associé à chaque cible.
- **QuantReports** (dans l'archive) : les résultats MassHunter du labo, qui servent de point de comparaison.

## 4. Le code

```text
src/
  agilent_ms.py        lecture des data.ms Agilent et export d'un CSV par ion (time_min, intensity)
  peaks.py             un pic : bruit, apex, bornes par double descente, aire
  reference.py         méthode (metadata.xlsx) et référence GAM-6 : TR et ratio de chaque composé
  measure.py           application de la référence à un échantillon (contrôles TR/ratio, réponse ISTD)
  run_reference.py     étape 1 : référence GAM-6
  run_batch.py         étape 2 : application de la référence aux 8 GAM et aux SF
  run_calibration.py   étapes 3-4 : calibration quadratique et contrôle SF
app.py                 application Streamlit d'analyse visuelle (lit outputs/)
tests/test_synthetic.py  validation sur chromatogrammes synthétiques à vérité connue
```

Pour lancer le pipeline :

```bash
python -m venv .venv && source .venv/bin/activate
python -m pip install -r requirements.txt
python src/run_reference.py data/raw
python src/run_batch.py data/raw
python src/run_calibration.py
streamlit run app.py
python tests/test_synthetic.py
```

## 5. Décisions prises

Chaque décision indique son statut : **validée** (vérifiée sur les données ou par le cours), **provisoire** (choix technique raisonnable, à confirmer) ou **à confirmer** (dépend du prof ou du cahier des charges).

| # | Décision | Justification | Statut |
| --- | --- | --- | --- |
| D1 | Extraire les ions directement des `data.ms` | Les CSV par ion n'existent pas dans l'archive. Lecteur vérifié identique à la bibliothèque `rainbow-api` sur des fichiers Agilent publics. | Validée |
| D2 | Un ion absent d'un scan SIM n'est **pas** mis à zéro | Il n'a pas été mesuré, ce n'est pas une intensité nulle | Validée |
| D3 | Bornes d'intégration par double descente depuis l'apex, avec un seuil en multiples du bruit estimé sur chaque trace | Méthode demandée en cours ; aucun temps codé en dur | Validée |
| D4 | Arrêt de la descente au retour à la ligne de base, ou à la vallée si un pic voisin remonte (séparation verticale des coélutions) | Gère les ions partagés (Benzo b/k, Phénanthrène/Anthracène…) | Provisoire, à comparer aux bornes MassHunter |
| D5 | Ligne de base horizontale au minimum des deux bornes | Convention de l'application montrée par le prof | Provisoire |
| D6 | Aires calculées sur le signal brut (trapèzes) ; lissage utilisé seulement pour localiser l'apex et les bornes | Ne pas déformer l'aire | Validée |
| D7 | Le qualifiant est intégré sur les bornes du quantifiant, resserrées aux vallées du qualifiant | Avec des bornes indépendantes, le ratio est sous-estimé quand le signal est bruité (jusqu'à −14 % en simulation) | Provisoire. L'appli du prof semble utiliser des bornes indépendantes. |
| D8 | Affectation des pics par ordre d'élution, avec contrôle croisé entre m/z | Seul invariant fiable. Ce contrôle a corrigé le Chrysène-D12, confondu avec le surrogat du m/z 240 avant l'ajout des surrogats. | Validée |
| D9 | Calibration quadratique pondérée en 1/x | Méthode affichée dans les rapports MassHunter | Validée par comparaison (voir 6) |
| D10 | Inversion : plus petite racine positive dans le domaine de la gamme ; hors domaine signalé, jamais extrapolé en silence | Cohérence physique | Validée |
| D11 | Tolérance SF fixée à ±20 % | Valeur de travail pour faire tourner le calcul | **À confirmer avec le cahier des charges** |
| D12 | Nominal SF = 1 ppm, facteur de dilution = 1 pour GAM et SF | Cours et inventaire du 5 octobre | À confirmer |

## 6. Résultats au 9 octobre 2026 (batch 20251103)

- **Référence** : 23 composés sur 23, chacun avec un TR et un ratio de référence. Les TR suivent exactement l'ordre de `metadata.xlsx`.
- **Gamme** : 181 mesures conformes sur 184. Les 3 exceptions sont le Fluorène et l'Anthracène en GAM-1 et GAM-2, à +0,22 à +0,25 % de TR, avec un ratio correct.
- **Constat sur le TR** :
  - le TR baisse régulièrement quand la concentration augmente, pour presque tous les composés (environ +0,25 % à 0,025 ppm, environ −0,1 % à 5 ppm) ;
  - le Pérylène-D12, qui ne coélue avec rien, ne bouge pas : c'est un effet de la concentration, pas de la colonne ;
  - c'est cohérent avec la remarque du cours (« le temps de passage dépend de la charge ») ;
  - l'exception est le Benzo(k)fluoranthène aux deux plus bas niveaux, où l'écart s'inverse (probable effet de la coélution avec le Benzo(b), à vérifier visuellement).
- **Calibration** : R² d'au moins 0,9997 pour les 16 HAP. Les réponses sont à ±5 % de celles de MassHunter et les coefficients b sont proches (Naphtalène 1,169 contre 1,195 ; Chrysène 1,261 contre 1,278). L'intégration est donc validée par le logiciel du labo.
- **SF** : la concentration retrouvée est de 0,88 à 0,99 ppm pour un nominal de 1 ppm, soit environ −9 % sur tous les HAP. Les 64 résultats sont PASS à ±20 %, mais cet écart systématique est à expliquer.
- **Étalons internes** : leurs aires restent dans ±7 % sur la gamme. Elles baissent sur les SF (0,85 à 0,95 de GAM-6), signe d'une perte de sensibilité au fil de la séquence, que la normalisation corrige.

## 7. Questions ouvertes (à poser au prof)

1. La tolérance de ±0,2 % sur le TR s'applique-t-elle strictement aux bas niveaux de gamme, vu la dérive liée à la concentration ?
2. Quelle tolérance d'acceptation pour les SF, et d'où vient l'écart systématique d'environ −9 % (nominal exact, préparation) ?
3. Le ratio qualifiant/quantifiant doit-il être intégré sur des bornes communes ou indépendantes (D7) ?
4. Quel est l'effet attendu de `double_peak` (Dibenz(a,h)anthracène, qui présente un épaulement, aujourd'hui inclus dans l'aire) et de `median_despike` (Indéno) ?
5. Les bornes de fin de pic sont larges sur les pics très intenses (jusqu'à 0,5 min de traînée). Faut-il les aligner sur MassHunter ? La comparaison bornes à bornes avec les QuantReports reste à faire.
6. Peut-on versionner `outputs/` (résultats calculés, quelques Mo, sans données brutes) pour déployer l'application sur Streamlit Cloud ?

## 8. Travailler à plusieurs (groupe et assistants)

- GitHub et ce fichier sont la mémoire commune. Une décision qui n'est écrite qu'en conversation n'existe pas.
- Avant de reprendre un chantier : lire ce fichier, `git log`, puis le code concerné.
- Toute décision structurante est ajoutée au tableau de la section 5, avec sa justification et son statut. Une décision existante n'est jamais supprimée : on la marque « remplacée par Dx ».
- Une contradiction entre une information orale, une ancienne analyse et le cahier des charges ou les données est signalée, pas tranchée en silence.
- Le travail passe par une branche et une PR relue par au moins un membre (voir `CONTRIBUTING.md`). Les audits d'assistant se font sur la PR.
- Les données brutes ne sont jamais versionnées.
