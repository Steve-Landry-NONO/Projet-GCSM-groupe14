# Analyse du rejeu GAM-6 — 9 octobre 2026

Sources : archives partage_rejeu_corrige.tar.gz et ions_gam6.tar.gz transmises
par le groupe ; classeur metadata(1).xlsx déjà identifié par son empreinte.
Les traces comprennent 33 ions. Aucune donnée brute n’est ajoutée au dépôt.

## Résultats du rejeu strict

- 276 lignes : 23 composés sur huit GAM et quatre SF.
- 192 mesures cibles : 128 GAM et 64 SF.
- 189/192 contrôles TR réussis et 192/192 contrôles de ratio réussis.
- Aucune réponse stricte : 16 calibrations impossibles et 64 SF sans calibration.
- Une alerte de vallée apparaît sur 221 lignes, y compris par propagation de la référence.

## Ce que montrent les signaux de GAM-6

Une vallée détectée n’implique pas systématiquement un chevauchement important.
Avec un fond descriptif estimé au dixième percentile de la trace, l’intensité
à la vallée des benzo(b)/benzo(k)fluoranthènes représente environ 47 et 52 %
de leurs hauteurs respectives. Plusieurs autres vallées sont beaucoup plus
basses : environ 0,8–1 % pour Phénanthrène/Anthracène et 1,5–1,9 % pour
Benzo(a)anthracène/Chrysène. Ce sont des indicateurs de forme, pas des erreurs
d’aire ni des seuils d’acceptation.

L’Acénaphtène-D10 atteint la fin de sa fenêtre SIM avec un signal résiduel
d’environ 0,03 % de la hauteur selon ce même indicateur. Cela ne prouve ni
une troncature importante ni une aire complète : la partie hors acquisition
est inconnue. Le message est donc « borne au bord de la fenêtre SIM : aire à examiner ».

## Correction de méthode logicielle

Le blocage automatique de toute vallée était trop général pour explorer
les résultats. On distingue désormais le calcul strict et le calcul exploratoire.

- Le calcul strict et ses exclusions restent disponibles.
- Une réponse exploratoire exige les contrôles TR/ratio de la cible et de son
  ISTD, des aires positives et finies et aucune attribution d’apex en double.
- Les alertes de référence, de vallée, de fenêtre ou de cas particulier restent
  visibles ; leur présence n’est pas une validation de l’intégration.
- La calibration exploratoire exige quatre concentrations distinctes et une
  courbe croissante, avec contrôle des deux bornes lors de l’inversion.
- Les SF exploratoires ne sont jamais PASS, quel que soit leur écart.

## Vérifications

19 tests de régression réussis, dont deux nouveaux vérifiant la séparation
des réponses stricte/exploratoire et l’absence de PASS SF en mode exploratoire.
Rejeu local des seules traces GAM-6 reçues : zéro réponse stricte, 16 réponses
exploratoires cibles. Les huit GAM et quatre SF n’ont pas encore été rejoués
avec cette nouvelle version, car leurs traces n’ont pas été transmises.

## Suite

Rejouer le batch dans un dossier séparé avec --exploratory. Examiner les points
exclus, les courbes, les résidus et les écarts SF. Comparer les bornes à
MassHunter avant toute acceptation quantitative. Les coélutions élevées et
les doubles pics demandent toujours une décision du professeur.
