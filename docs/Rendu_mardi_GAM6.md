# Rendu mardi — Références GAM-6

## Objectif

À partir du 6e point de gamme du batch `20251103` (`20251103-GAM-25-385-6.D`), déterminer pour les 16 HAP cibles :

- le temps de rétention de référence, pris à l'apex du pic quantifiant ;
- les fenêtres d'intégration des ions quantifiant et qualifiant ;
- le ratio de référence `aire_qualifiant / aire_quantifiant`.

## Méthode

1. Extraction des traces ioniques depuis `data.ms` vers des séries `time_min / intensity`.
2. Détection des pics par proéminence par rapport au bruit.
3. Pour les m/z partagés, affectation des pics selon l'ordre des composés défini dans `metadata.xlsx`.
4. Recherche des bornes de chaque pic par descente depuis l'apex jusqu'au retour à la baseline ou à une vallée entre deux pics.
5. Intégration sur le signal brut, avec une baseline horizontale au plus bas des deux bords.
6. Calcul indépendant des fenêtres quantifiant et qualifiant, puis calcul du ratio des aires.

Aucune tolérance métier (±0,2 %, ±23 %, etc.) n'est appliquée dans ce rendu : ces valeurs doivent être confirmées par une source métier avant d'être utilisées comme critères de conformité.

## Résultats historiques du premier calcul

Ces valeurs sont conservées pour tracer le travail initial. Elles n’ont pas été recalculées lors de la correction du 9 octobre. Elles ne constituent pas une validation de la quantification. Le contrôle « OK » signifie ici absence d’alerte rapportée dans le premier rendu ; les coélutions et doubles pics restent à examiner.

## Tableau du premier rendu

| Composé | m/z quant | m/z qual | RT réf. (min) | Ratio qual/quanti | Contrôle |
|---|---:|---:|---:|---:|---|
| Naphtalene | 128 | 127 | 4.8358 | 0.1278 | OK |
| Acenaphtylene | 152 | 153 | 5.8782 | 0.1301 | affectation à vérifier |
| Acenaphtene | 153 | 154 | 5.9874 | 0.9341 | borne à vérifier |
| Fluorene | 166 | 165 | 6.3431 | 0.8975 | OK |
| Phenanthrene | 178 | 176 | 7.0396 | 0.1823 | OK |
| Anthracene | 178 | 176 | 7.0964 | 0.1745 | OK |
| Fluoranthene | 202 | 101 | 8.4647 | 0.1319 | OK |
| Pyrene | 202 | 101 | 8.8178 | 0.1563 | OK |
| Benzo(a)anthracene | 228 | 226 | 10.8392 | 0.2536 | OK |
| Chrysene | 228 | 226 | 10.9141 | 0.2762 | OK |
| Benzo(b)fluoranthene | 252 | 253 | 12.6556 | 0.2151 | coélution proche |
| Benzo(k)fluoranthene | 252 | 253 | 12.6863 | 0.2138 | coélution proche |
| Benzo(a)pyrene | 252 | 253 | 13.1762 | 0.2156 | OK |
| Dibenz(a,h)anthracene | 278 | 279 | 14.6815 | 0.2250 | double_peak à confirmer |
| Indeno(1,2,3-cd)pyrene | 276 | 277 | 14.7300 | 0.2338 | affectation à vérifier |
| Benzo(g,h,l)perylene | 276 | 277 | 15.0812 | 0.2350 | affectation à vérifier |

## Points à vérifier

- `Benzo(b)fluoranthene` et `Benzo(k)fluoranthene` sont très proches sur m/z 252/253 : la séparation se fait à la vallée entre les deux pics.
- `Dibenz(a,h)anthracene` est marqué `double_peak` dans les metadata : le pic retenu ici est le pic principal ; ce cas doit être confronté à la règle métier attendue.
- Les avertissements d'affectation indiquent qu'un contrôle visuel reste recommandé avant de figer définitivement les références.


## Reprendre ce travail

Ce module est le prototype historique, avec des bornes qualifiantes indépendantes.
Le pipeline de la PR #2 utilise une autre convention configurable : les résultats
ne doivent pas être mélangés sans conserver le nom et la version de la méthode.
Après fusion des deux PR, le point d'entrée courant sera `src/run_reference.py`.

Installer les dépendances avec `python -m pip install -r requirements-gam6.txt`.
Vérifier les cas élémentaires avec :

```bash
python -m unittest discover -s tests -p test_gam6_reference.py -v
```

Les corrections conservent une ligne et une alerte si le quantifiant manque,
si le qualifiant manque ou si sa fenêtre ne recouvre pas celle du quantifiant.
Le calcul d'aire fonctionne avec NumPy 1.x et 2.x. Le lecteur Agilent est identique
sur les deux branches pour préparer leur fusion.

Le [plan global](AUDIT_ET_PLAN_GLOBAL.md) décrit les limites et la suite du projet.
