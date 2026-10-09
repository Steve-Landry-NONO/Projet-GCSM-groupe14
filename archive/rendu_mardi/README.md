# Archive — premier rendu GAM-6 (rendu du mardi)

Ce dossier conserve le premier calcul des références GAM-6, rendu avant le pipeline actuel. Il est gardé pour la traçabilité du travail. **Ce n'est pas un point d'entrée du projet.** Le calcul courant passe par `src/run_reference.py` (voir le README à la racine).

## Contenu et provenance

| Fichier | Emplacement d'origine | Historique (branche `rendu-mardi-gam6`, PR #1) |
| --- | --- | --- |
| `gam6_reference.py` | `src/gam6_reference.py` | Créé le 7 octobre 2026 (450d29e), sécurisé le 9 octobre (8fc20db) |
| `gam6_reference_20251103.csv` | `results/gam6_reference_20251103.csv` | Produit et versionné le 7 octobre 2026 (a483117) |
| `Rendu_mardi_GAM6.md` | `docs/Rendu_mardi_GAM6.md` | Écrit le 7 octobre (658dffd), complété le 9 octobre (de7d931) |
| `test_gam6_reference.py` | `tests/test_gam6_reference.py` | Ajouté le 9 octobre (8fc20db) |
| `requirements-gam6.txt` | `requirements-gam6.txt` | Ajouté le 9 octobre (8fc20db) |

Les fichiers ont été déplacés avec `git mv` : leur historique reste consultable avec `git log --follow archive/rendu_mardi/<fichier>`.

Le lecteur Agilent de ce prototype n'est pas archivé : `src/agilent_ms.py` est commun au prototype et au pipeline, et identique sur les deux branches.

## Ce que représentent ces résultats

- **Données** : batch `20251103`, acquisition `20251103-GAM-25-385-6.D`, 16 HAP cibles. Les étalons internes et les surrogats ne sont pas inclus.
- **Méthode du prototype** : bornes par descente depuis l'apex, ligne de base horizontale au plus bas des deux bornes, et **fenêtres qualifiantes indépendantes** du quantifiant. Le pipeline courant intègre le qualifiant sur les bornes du quantifiant (décision D7 de `CONTEXTE_PROJET.md`). Les ratios des deux méthodes ne sont donc pas directement comparables.
- **Statut** : résultats historiques, produits le 7 octobre avec la première version du module. Ils n'ont pas été recalculés après les corrections du 9 octobre. Ils ne valident aucune quantification.

## Relancer les tests du prototype

Depuis la racine du dépôt :

```bash
python -m pip install -r archive/rendu_mardi/requirements-gam6.txt
python -m unittest discover -s archive/rendu_mardi -p test_gam6_reference.py -v
```
