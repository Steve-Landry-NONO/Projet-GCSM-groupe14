# Projet GCSM — Groupe 14

## Analyse de données GC-MS pour la quantification de molécules toxiques

Ce dépôt rassemble le travail du groupe 14 sur l'exploitation de données de chromatographie en phase gazeuse couplée à la spectrométrie de masse (GC-MS). Le projet vise à déterminer les concentrations de molécules toxiques connues, ici 16 hydrocarbures aromatiques polycycliques (HAP), à partir des signaux instrumentaux, d'une gamme d'étalonnage et de solutions de contrôle.

Le dépôt est volontairement évolutif. Il contient le cadrage, un notebook d'exploration et, depuis le 9 octobre 2026, un pipeline Python qui recalcule la référence GAM-6, les mesures de la gamme, les calibrations et le contrôle des SF à partir des signaux ioniques bruts, ainsi qu'une application Streamlit pour analyser les résultats.

Pour reprendre le projet, commencer par [CONTEXTE_PROJET.md](CONTEXTE_PROJET.md) : contexte métier, décisions prises, résultats et questions ouvertes.

## Membres du groupe

- Steve Landry KOUOKAM — chef de groupe — `@Steve-Landry-NONO`
- Ludovic TUEKAM — `@ludovictuekam9-hue`
- Radia GHILAS — `@radiaghilas`
- Harald MAFORIKAN — `@harald8`
- Ismaila DIEYE — `@ismailadieye-iage2`

Encadrement et relecture : `@HorHakim` et `@Septentrion`.

## Question centrale

Quelles informations devons-nous extraire des données GC-MS et sous quelle forme devons-nous les organiser pour identifier les molécules recherchées et calculer leurs concentrations de manière traçable ?

Notre réponse actuelle est la suivante : un nom de composé ou un chromatogramme global ne suffit pas. Pour chaque molécule, il faut réunir les ions quantifiant et qualifiant, le temps de rétention, les aires des pics, l'étalon interne associé, les concentrations connues des étalons, le facteur de dilution et les règles de validation. La concentration n'est exploitable qu'après calibration et vérification sur une solution de contrôle.

## Livrables actuels

| Livrable | Rôle |
| --- | --- |
| [Rapport PDF](docs/Rendu_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.pdf) | Version de référence à remettre, avec mise en page fixe, police Calibri et figures en noir et blanc. |
| [Rapport Markdown](docs/Rendu_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.md) | Version lisible directement sur GitHub, facile à corriger et à faire évoluer. |
| [Notebook Jupyter](notebooks/Notebook_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.ipynb) | Exploration reproductible, commentaires simples et visualisations. |
| [Notebook HTML](notebooks/Notebook_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.html) | Aperçu autonome du notebook déjà exécuté. |
| [Figures](https://github.com/Steve-Landry-NONO/Projet-GCSM-groupe14/tree/main/assets/figures) | Images en noir et blanc utilisées dans le rapport et dans le notebook. |
| [Pipeline `src/`](src) | Extraction des ions depuis les `.D`, référence GAM-6, mesures GAM/SF, calibration quadratique et contrôle SF. |
| [Application `app.py`](app.py) | Analyse visuelle interactive : référence, gamme, calibration, contrôle SF. |
| [CONTEXTE_PROJET.md](CONTEXTE_PROJET.md) | Contexte, décisions, résultats et questions ouvertes, à lire avant toute reprise. |

Le PDF conserve exactement la présentation du rendu. Le Markdown conserve le contenu scientifique, les tableaux, les formules et les liens vers les figures, mais GitHub applique sa propre police et sa propre mise en page. Les deux formats sont donc complémentaires.

## Ce que nous avons compris des données

L'archive explorée contient 46 acquisitions instrumentales :

| Type d'acquisition | Nombre | Utilité comprise |
| --- | ---: | --- |
| GAM | 8 | Construire une courbe d'étalonnage sur huit niveaux : 0,025 à 5 ppm. |
| SF | 4 | Contrôler la justesse de la calibration sur une solution connue. |
| Hexane | 6 | Contrôler le solvant ou un blanc ; le rôle exact reste à confirmer. |
| BLPC | 28 | Appliquer ensuite la méthode à des échantillons de terrain. |

La chaîne de traitement retenue est :

1. construire les repères de référence à partir de GAM-6 ;
2. repérer et intégrer les pics dans tous les niveaux GAM ;
3. normaliser l'aire de chaque composé par son étalon interne ;
4. ajuster une calibration quadratique, avec la pondération à confirmer ou à reproduire ;
5. vérifier cette calibration sur les solutions SF ;
6. appliquer la méthode validée aux solutions inconnues et corriger par le facteur de dilution.

## Structure de données proposée

Nous proposons de garder les fichiers instrumentaux natifs sans modification et de construire des tables liées par `batch_id`, `sample_id` et `compound_id`.

| Table | Contenu principal |
| --- | --- |
| `samples.csv` | Identité de l'acquisition, type, niveau GAM, concentration connue, date et dilution. |
| `compounds.csv` | Molécules, ions, étalon interne, ordre des pics et cas particuliers. |
| `signals.parquet` | Temps et intensité pour chaque ion mesuré. |
| `references.csv` | Temps de rétention et fenêtres d'intégration de référence. |
| `peaks.parquet` | Aires, ratios, bornes d'intégration et anomalies. |
| `calibrations.csv` | Coefficients, domaine, pondération, niveaux utilisés et diagnostics. |
| `results.csv` | Concentrations avant et après dilution, unité, statut et justification. |

## Organisation du dépôt

```text
Projet-GCSM-groupe14/
├── README.md
├── CONTEXTE_PROJET.md
├── CONTRIBUTING.md
├── requirements.txt
├── app.py                      application Streamlit
├── src/                        pipeline de traitement
│   ├── agilent_ms.py           lecture des .D
│   ├── peaks.py                détection et intégration d'un pic
│   ├── reference.py            référence GAM-6
│   ├── measure.py              contrôle des GAM et SF
│   ├── alerts.py               alertes bloquantes ou informatives
│   ├── run_reference.py
│   ├── run_batch.py
│   └── run_calibration.py
├── tests/                      tests unitaires et caractérisation synthétique
├── archive/rendu_mardi/        prototype du premier rendu GAM-6, historique (apporté par la PR #1)
├── data/                       (non versionné) data/raw/<batch>/*.D et data/metadata.xlsx
├── outputs/                    (non versionné) résultats calculés
├── docs/
│   ├── Rendu_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.md
│   └── Rendu_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.pdf
├── notebooks/
│   ├── Notebook_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.ipynb
│   └── Notebook_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.html
└── assets/
    └── figures/
```

## Lancer le pipeline et l'application

Placer les dossiers `.D` du batch dans `data/raw/` et `metadata.xlsx` dans `data/`, puis :

```bash
python -m venv .venv
source .venv/bin/activate          # Windows : .venv\Scripts\activate.bat
python -m pip install -r requirements.txt

python src/run_reference.py data/raw/20251103     # 1. référence GAM-6 : TR et ratio de référence
python src/run_batch.py data/raw/20251103         # 2. mesures des 8 GAM et des SF par rapport à la référence
python src/run_calibration.py            # 3-4. calibration quadratique (1/x) et contrôle SF
streamlit run app.py                     # analyse visuelle
python tests/test_synthetic.py           # validation sur signaux synthétiques
```

Les résultats sont écrits dans `outputs/<batch>/` : `references_gam6.csv`, `peaks.csv`, `conformite.csv`, `calibrations.csv`, `calibration_points.csv`, `sf_results.csv`, ainsi que des figures de contrôle.

## Reproduire le notebook

Créer un environnement Python, installer les dépendances puis ouvrir le notebook :

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
jupyter notebook notebooks/Notebook_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.ipynb
```

Le notebook contient les données de synthèse nécessaires à ses graphiques. Les fichiers GC-MS bruts ne sont pas ajoutés à ce dépôt, car ils sont volumineux et doivent rester associés à leur provenance instrumentale.

## État et prochaines étapes

- [x] Explorer l'archive et inventorier les acquisitions.
- [x] Identifier les 16 HAP, leurs ions et leurs étalons internes.
- [x] Décrire les données indispensables et proposer un modèle exploitable.
- [x] Produire un PDF, un Markdown et un notebook commenté.
- [x] Extraire les signaux séparés par ion depuis les dossiers `.D`.
- [x] Construire les références à partir de GAM-6 (23 composés).
- [x] Implémenter la mesure de la gamme GAM ; résultats à recalculer après les garde-fous.
- [x] Implémenter les calibrations ; leur validité dépend désormais des points sans alerte.
- [ ] Valider les SF : calcul fait, tolérance et écart systématique d'environ −9 % à confirmer.
- [ ] Confirmer les règles d'intégration et les tolérances métier (voir les questions ouvertes de CONTEXTE_PROJET.md).
- [ ] Quantifier les échantillons inconnus (BLPC) et documenter les incertitudes.

## Points à confirmer

- concentration nominale exacte des solutions SF ;
- tolérances d'acceptation sur le temps de rétention, le ratio ionique et la concentration ;
- méthode officielle d'intégration et de ligne de base ;
- interprétation de `double_peak`, `median_despike` et de la coélution ;
- ~~moyen retenu pour convertir les données natives `.D` en séries `temps/intensité` par ion~~ : lecture directe des `data.ms` (`src/agilent_ms.py`).

## Règles de contribution

Les règles simples de travail, les noms de branches et la convention de messages de commit sont décrits dans [CONTRIBUTING.md](CONTRIBUTING.md). Les fichiers sources modifiables doivent être mis à jour avant les exports PDF ou HTML correspondants.


## Audit et plan commun des branches

Le [plan global](docs/AUDIT_ET_PLAN_GLOBAL.md) compare `main`, le premier rendu
GAM-6 (PR #1) et le pipeline (PR #2). Les résultats historiques restent consultables,
mais doivent être régénérés avant une nouvelle conclusion scientifique.

Le traitement refuse plusieurs batchs à la fois. Une alerte sur le pic, la référence
ou l’ISTD bloque l’usage automatique de la réponse en calibration. Les points exclus
restent dans `calibration_points.csv`, avec leur motif. Le statut `ok` d’une courbe
indique ici un ajustement calculable et croissant, pas une validation complète du laboratoire.
Les coélutions, doubles pics et filtres de pics parasites demandent encore une validation métier.

```bash
python -m unittest discover -s tests -p "test_*.py" -v
python tests/test_synthetic.py
```

L’application affiche les graphiques en noir et blanc, avec Calibri si la police est
installée (Arial en repli). Le succès du test synthétique concerne les aires isolées
et le signalement des limites : il ne valide pas les aires coéluées.

## Calcul exploratoire après le rejeu

Voir [l’analyse de GAM-6](docs/ANALYSE_REJEU_GAM6.md). Une vallée ou une fin de
fenêtre reste une alerte, sans prouver à elle seule une erreur d’aire.
`response_exploratory` permet un diagnostic si la cible et son ISTD passent
les contrôles TR/ratio et ont des aires positives et finies. Les alertes restent
visibles ; la réponse stricte reste bloquée.

`run_calibration.py --exploratory` utilise cette réponse uniquement pour comparer
les courbes. Les calibrations sont marquées `exploratoire` et toutes les SF
calculées `NON VALIDÉ (exploratoire)`, même avec un faible écart. Utiliser un
dossier de sortie distinct des résultats stricts. Les seuils ne sont pas élargis.

## Alertes bloquantes et informatives

Une vallée entre deux pics n'est plus bloquante par principe. Le rapport
vallée/hauteur `v = (L_borne − F)/(H − F)` est mesuré sur chaque pic ; la
vallée devient informative si `v <= 10 %` (`--valley-max`), et un bord de
fenêtre SIM si le signal résiduel est `<= 1 %` de la hauteur (`--edge-max`).
Toutes les autres alertes restent bloquantes : attribution ambiguë, pic
réaffecté, qualifiant douteux, `double_peak`, `median_despike`. Ces seuils
sont expérimentaux, à confirmer avec l'encadrement : ils ne prouvent pas la
justesse des aires. Méthode, tests et limites :
[docs/VALLEE_HAUTEUR.md](docs/VALLEE_HAUTEUR.md).

```bash
python tests/caracteriser_vallee.py      # tableau de caractérisation (cas synthétiques)
```
