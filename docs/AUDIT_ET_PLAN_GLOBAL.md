# Audit et plan global — GC-MS, groupe 14

Date : 9 octobre 2026.

## Pourquoi ce document

Nous voulons comprendre les deux travaux en cours, les réunir sans perdre le premier rendu et continuer à partir de résultats dont les limites sont explicites. Une aire, un ratio et une concentration calculée ne représentent pas le même contrôle.

## État comparé avant correction

| Branche | Version auditée | Contenu | Rôle après réunion |
| --- | --- | --- | --- |
| main | 1e8dfba | Rapport PDF/Markdown, notebook et figures du cadrage | Conserver les livrables initiaux |
| rendu-mardi-gam6 — PR #1 | 450d29e | Premier lecteur Agilent, référence GAM-6, CSV de 16 cibles et note du rendu | Archive du premier calcul et de sa méthode |
| analysis/reference-gam6 — PR #2 | 5b3825a | Référence de 23 composés, GAM/SF, normalisation ISTD, calibration, application | Pipeline courant après vérification et fusion |

Les deux PR partaient du même main. Le lecteur src/agilent_ms.py était ajouté dans les deux branches avec des différences de présentation et un lancement en ligne de commande dans la PR #2. Son contenu est désormais harmonisé. Les autres modules de référence ont des noms et des conventions différents : ils restent distingués.

Le notebook présent sur main est une exploration avec données embarquées, TIC et diagnostics publiés. Il ne calcule pas lui-même les concentrations à partir des ions. Le README de main décrit correctement ce premier stade mais ne représente pas encore les avancées des branches. Son évolution est préparée dans la PR #2. Le notebook local lu pendant l'audit a été comparé à l'empreinte GitHub : e23dfda2b0decbb741dec9089fb866f022071ddb.

## Ce que nous calculons

1. Une série par ion : temps en minutes et intensité mesurée.
2. Une référence GAM-6 : temps d'apex, bornes et ratio qualifier/quantifier.
3. Dans chaque acquisition : une identification et une aire avec leurs alertes.
4. Une réponse : aire du quantifier de la cible / aire du quantifier de son ISTD.
5. Une courbe : réponse = a C² + b C + c.
6. Une concentration obtenue par inversion, avec contrôle du domaine.
7. Un contrôle SF utilisant le statut du pic, la calibration et l'écart au nominal.

Les ratios d'identification et de normalisation sont différents. Un bon ratio d'identification ne suffit pas à prouver une bonne aire.

## Corrections de la PR #1

- Signaler un quantifiant ou un qualifiant absent au lieu de provoquer une erreur.
- Signaler une fenêtre qualifiante sans recouvrement exploitable.
- Accepter une liste vide de composés.
- Rendre l'intégrale compatible avec NumPy 1.x et 2.x.
- Harmoniser le lecteur Agilent avec la PR #2.
- Ajouter cinq tests élémentaires et requirements-gam6.txt.
- Conserver les nombres du CSV comme résultats historiques non régénérés.

Ce prototype conserve des bornes qualifiantes indépendantes. Il ne doit pas devenir un deuxième pipeline de production concurrent.

## Corrections de la PR #2

- Refuser une source qui contient plusieurs parents batchs.
- Garder une ligne pour les composés dont la référence manque.
- Tenir compte de use et applicable_to dans les mesures ; le champ Inutilisé et les composés non applicables restent explicitement signalés.
- Propager les alertes de référence et d'intégration.
- Signaler les vallées entre pics, doubles pics et corrections median_despike encore non définies.
- Bloquer la normalisation si la cible ou l'ISTD est invalide ; exiger des aires positives et finies.
- Signaler un même apex attribué à plusieurs composés du même m/z.
- Exclure les points douteux de la calibration tout en conservant le motif dans calibration_points.csv.
- Exiger au moins quatre concentrations distinctes exploitables.
- Refuser une courbe non croissante ; sélectionner la racine à pente positive.
- Contrôler la borne basse et la borne haute ; toute extrapolation est marquée hors gamme.
- Ne jamais transformer un pic douteux en PASS SF à partir du seul écart de concentration.
- Réutiliser le même calcul de verdict dans le script et dans l'application.
- Retirer le cache fondé uniquement sur le chemin des fichiers pour éviter les anciens résultats après recalcul.
- Passer les graphiques de l'application en niveaux de gris, sans emoji ; demander Calibri avec Arial en repli.
- Gérer des résultats vides ou insuffisants sans afficher une ancienne figure de calibration.

Ces garde-fous sont conservateurs : ils bloquent certains résultats qui demandaient déjà une relecture. Une procédure de validation manuelle documentée reste à concevoir.

## Vérifications effectuées

- PR #1 : cinq tests unittest réussis.
- PR #2 : dix-sept tests unittest réussis, y compris deux exécutions de scripts sur des dossiers temporaires.
- Simulation : aires des pics isolés retrouvées à moins de 5 % sur les exemples testés ; contrôles de répétabilité du TR et du ratio sur 200 tirages.
- Coélutions simulées : erreurs d'aire allant jusqu'à environ −35,8 % et +37,0 %. Elles restent présentes et entraînent une alerte. Le test ne prétend plus les valider.
- Compilation Python des scripts et de l'application.

Les acquisitions natives du laboratoire n'ont pas été rejouées lors de cet audit. Les valeurs historiques 181/184, R² >= 0,9997 et 64/64 PASS ne décrivent donc pas les résultats après correction. La compilation de l'application ne remplace pas un test de son interface dans un environnement Streamlit équipé des sorties réelles.

## Données et documents à réconcilier

Le metadata.xlsx joint contient 23 composés : 16 cibles, 4 ISTD et 3 surrogates.

Empreinte SHA-256 de la pièce jointe metadata(1).xlsx :
dc0afab6dd6860450db8aedc684fa4e20c74596f0d563249edc2a5600250beb1

Le classeur reste externe au dépôt. Il faut utiliser exactement cette version pour la prochaine exécution, puis consigner son empreinte avec les résultats.

Le nouvel Excalidraw ne figure pas dans les fichiers modifiés des PR auditées ; sa version courante reste à récupérer. Les notes parlent de 29 séries et le contexte de 33 ions : réconcilier les ions réellement exportés et les ions utiles à la méthode.

## Plan de continuation

| Étape | Travail | Critère de fin |
| --- | --- | --- |
| 1 — Sécurisation | Corriger les erreurs logicielles et tracer les exclusions | Tests de régression réussis ; corrections relues |
| 2 — Validation métier | Confirmer seuils TR/ratio/SF, nominal SF, volumes, doubles pics, despiking, bornes du qualifier | Décisions écrites avec la source ou la validation du professeur |
| 3 — Rejeu réel | Refaire GAM-6, GAM et SF avec le nouveau classeur | Tables et figures de contrôle, exclusions expliquées |
| 4 — Intégration | Comparer les bornes et les aires à MassHunter, particulièrement pour les coélutions | Écarts mesurés ; règle acceptée ou méthode de séparation revue |
| 5 — Calibration | Examiner les points exclus, résidus, bas niveaux, monotonie et domaine | Critères d'acceptation définis et appliqués |
| 6 — SF | Vérifier nominal et préparation, puis expliquer le biais commun | Contrôles acceptés selon des règles confirmées |
| 7 — BLPC | Lire la dilution, contrôler les blancs et surrogates, calculer les inconnus | Résultats traçables avec statut, unités et limites |
| 8 — Terre | Rassembler masse extraite, volume final, humidité éventuelle et dilution | Conversion solution vers mg/kg justifiée |
| 9 — Rendu | Actualiser le notebook et les exports PDF/HTML après validation des sources | Résultats et explications cohérents dans tous les formats |

## Préparer les merges

1. Faire relire les corrections par un membre du groupe, conformément à CONTRIBUTING.md.
2. Fusionner la PR #1 pour conserver le premier rendu.
3. Vérifier la PR #2 contre le nouveau main et exécuter aussi les tests du prototype.
4. Fusionner la PR #2 lorsque les garde-fous sont acceptés et les limites scientifiques clairement assumées.
5. Relire le README sur main et refaire un lancement dans un environnement propre.

Aucun merge n'est réalisé pendant cette correction. L'accord sur une fusion de code ne signifie pas que toutes les concentrations sont validées.

## Questions prioritaires au professeur

- Les seuils ±0,2 % et ±23 % s'appliquent-ils aussi aux points GAM les plus faibles ?
- Pour une coélution, faut-il couper à la vallée, ajuster plusieurs profils, ou reprendre les bornes MassHunter ?
- Les aires qualifiantes utilisent-elles des bornes propres ou celles du quantifiant ?
- Quelle règle exacte pour double_peak, median_despike et coelution_order ?
- Quel nominal SF, quelle dilution et quels critères de calibration/SF officiels ?
- Quelle unité finale est attendue, avec quels volumes et quelle masse de terre ?

## Repères chimiques

Prélever une partie d'une solution homogène ne la dilue pas. Le naphtalène-d8 contient huit deutériums à la place des huit hydrogènes du naphtalène. L'ionisation électronique provoque l'ionisation et la fragmentation par un faisceau d'électrons, pas par un second four à 2 500 °C.

Sources de vérification : [NIST — Naphthalene-D8](https://webbook.nist.gov/cgi/cbook.cgi?ID=1146-65-2) et [Agilent — GC/MS FAQs](https://www.agilent.com/en/product/gas-chromatography-mass-spectrometry-gc-ms/gcms-fundamentals/gcms-faqs).

