# Guide de contribution

## Principe général

Chaque modification doit être compréhensible par un autre membre du groupe. Un commit correspond à une idée principale : documentation, données, analyse, correction ou organisation du dépôt.

## Travail avec les branches

Pour une modification importante, créer une branche courte à partir de `main` :

- `docs/nom-court` pour le rapport ou le README ;
- `analysis/nom-court` pour le notebook et les calculs ;
- `data/nom-court` pour les schémas ou scripts de préparation ;
- `fix/nom-court` pour une correction.

## Messages de commit

Utiliser un verbe à l'infinitif ou un résumé clair après un préfixe :

```text
docs: préciser les données nécessaires à la quantification
analysis: ajouter les diagnostics de la gamme GAM
data: définir la table des signaux ioniques
fix: corriger la concentration nominale des SF
chore: réorganiser les livrables du groupe
```

Le corps du commit doit expliquer, si nécessaire, ce qui a été compris ou corrigé et pourquoi cela compte pour le calcul final.

## Vérifications avant une proposition de modification

1. Relire les noms des composés, les unités et les concentrations.
2. Exécuter toutes les cellules du notebook dans l'ordre.
3. Vérifier que les graphiques restent lisibles en noir et blanc.
4. Vérifier les liens du README et du rapport Markdown.
5. Mettre à jour le PDF ou le HTML si leur source a changé.
6. Ne pas ajouter de données brutes confidentielles ou très volumineuses sans accord du groupe.

## Relecture

Une modification scientifique doit être relue par au moins un autre membre du groupe. Les professeurs peuvent recevoir un accès en lecture ou en collaboration selon les règles choisies par le propriétaire du dépôt.

