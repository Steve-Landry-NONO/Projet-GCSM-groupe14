# Projet GCSM — Groupe 14

## Analyse de données GC-MS pour la quantification de molécules toxiques

Ce dépôt rassemble le travail du groupe 14 sur l'exploitation de données de chromatographie en phase gazeuse couplée à la spectrométrie de masse (GC-MS). Le projet vise à déterminer les concentrations de molécules toxiques connues, ici 16 hydrocarbures aromatiques polycycliques (HAP), à partir des signaux instrumentaux, d'une gamme d'étalonnage et de solutions de contrôle.

Le dépôt est volontairement évolutif. Il contient le cadrage actuel, les données jugées nécessaires, la structure de données proposée et un notebook de visualisation. À ce stade, il ne présente pas encore un recalcul complet des concentrations à partir des signaux ioniques bruts.

## Membres du groupe

- Steve Landry KOUOKAM — chef de groupe — `@Steve-Landry-NONO`
- Ludovic TUEKAM — `@ludovictuekam9-hue`
- Radia GHILAS — `@Radiaghilas`
- Harald MAFORIKAN — `@harald8`
- Ismaila DIEYE

Encadrement et relecture : `@Horhakim` et `@Septentrion`.

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
| [Figures](assets/figures/) | Images en noir et blanc utilisées dans le rapport et dans le notebook. |

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
├── CONTRIBUTING.md
├── requirements.txt
├── docs/
│   ├── Rendu_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.md
│   └── Rendu_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.pdf
├── notebooks/
│   ├── Notebook_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.ipynb
│   └── Notebook_KOUOKAM_TUEKAM_GHILAS_MAFORIKAN_DIEYE.html
└── assets/
    └── figures/
```

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
- [ ] Extraire ou obtenir les signaux séparés par ion depuis les dossiers `.D`.
- [ ] Définir précisément les règles d'intégration et les tolérances métier.
- [ ] Construire les références à partir de GAM-6.
- [ ] Recalculer les calibrations pour chaque composé.
- [ ] Valider les résultats avec les SF.
- [ ] Quantifier les échantillons inconnus et documenter les incertitudes.

## Points à confirmer

- concentration nominale exacte des solutions SF ;
- tolérances d'acceptation sur le temps de rétention, le ratio ionique et la concentration ;
- méthode officielle d'intégration et de ligne de base ;
- interprétation de `double_peak`, `median_despike` et de la coélution ;
- moyen retenu pour convertir les données natives `.D` en séries `temps/intensité` par ion.

## Règles de contribution

Les règles simples de travail, les noms de branches et la convention de messages de commit sont décrits dans [CONTRIBUTING.md](CONTRIBUTING.md). Les fichiers sources modifiables doivent être mis à jour avant les exports PDF ou HTML correspondants.
