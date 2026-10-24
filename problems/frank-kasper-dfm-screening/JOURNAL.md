# Journal de reprise — Frank–Kasper

## 2026-10-24 — 1. État initial et objectif

Reprise demandée par le propriétaire : viser une résolution complète, sans
sous-agents, et documenter progressivement l'avancement. Travail effectué par
Codex (GPT-6 ; effort non enregistré). Aucune vérification par un autre agent
n'est effectuée durant cette reprise.

Le dépôt était propre au début. Le dossier contient un moteur de minimisation,
un tableau de 124 structures et six maillages entiers accompagnés de certificats
de bornes supérieures. Il ne contient ni preuve de minimalité globale, ni
manuscrit LaTeX, ni rapport de vérification. Les anciens checkpoints et sources
ignorés par Git ne sont pas présents dans cette copie.

**Point logique central :** une borne supérieure pour σ, même inférieure aux
valeurs calculées pour tous les concurrents, ne prouve pas que leurs infima
sont supérieurs à celui de σ. Il faut des bornes inférieures sur les classes
entières de partitions admissibles, ou un contre-exemple certifié.

Étapes engagées :

1. Reprendre les définitions et contrôler les sources primaires disponibles.
2. Rejouer et auditer les six certificats, y compris la topologie et les boîtes.
3. Examiner les voies de preuve globale et les dégénérescences admissibles.
4. Consigner des résultats reproductibles, avec une séparation explicite entre
   preuves, calculs numériques et questions non résolues.

Le catalogue de Downloads/PDF a été inspecté : les fichiers présents concernent
les marches auto-évitantes, EDP, bifurcations et transport. Aucun titre
Frank–Kasper, Reddy ou Dutour n'a été trouvé dans cet inventaire.

## 2026-10-24 — 2. Sources et reproduction

- Les six certificats historiques passent de nouveau :
  `certificates/logs/recheck_2026-10-24.log`.
- Le contrôle Milestone 1, y compris les recalculs des mousses de Kelvin et
  Weaire–Phelan, termine avec zéro échec :
  `certificates/logs/milestone1_recheck_2026-10-24.log`.
- L'archive Zenodo a été restaurée ; son MD5 est bien
  `582371cee665c62e01a2f47a73dd2c9d`. Le contrôle historique du catalogue
  donne 121 PASS et zéro échec. Les trois énergies indicatives de ressort
  II:Y04–Y06 diffèrent légèrement lors de cette reproduction ; ce sont des
  optimisations flottantes, sans effet sur les six certificats d'énergie.
- L'utilisateur a fourni le PDF de l'article de 2012,
  DOI `10.1007/s11224-012-0005-3`. Introduction et tableaux 2–4 lus ; pages
  1103 et 1107 également contrôlées visuellement. Les faces courbes sont bien
  autorisées. Les séries N01–N26, E01–E05 et Y01–Y06 correspondent aux trois
  tableaux. L'exhaustivité annoncée est bornée à 21 cellules dans un domaine
  fondamental réduit, pas à toutes les structures FK.

Corrections logiques repérées :

- « Non-Voronoï avec les coordonnées fournies » ne signifie pas « impossible
  à réaliser comme diagramme de Voronoï ».
- Le nombre de réalisations de Voronoï dans le catalogue est 105, non 108 ;
  108 inclut les trois comparateurs supplémentaires P, M et BCC.
- Les valeurs des huit relaxations dégénérées sont des valeurs de départ,
  pas des minima. Aucune exclusion globale ne découle de leur hauteur.
- Le test historique du catalogue utilisait des tolérances flottantes pour
  l'immersion. Un vérificateur exact distinct est en préparation.

Un second programme de vérification des six maillages est en cours : contrôle
des sphères cellulaires, disques des faces et liens de sommets, identité du
complexe dual périodique, intégration rationnelle indépendante, puis transport
affine certifié vers les boîtes cristallographiques exactes.

## 2026-10-24 — 3. Nouveaux résultats certifiés

**Les 121 structures disposent désormais de certificats géométriques exacts.**
`scripts/certify_catalogue.py` vérifie sur une grille dyadique : déterminants
positifs, annulation orientée des faces, volume total exactement égal à celui
du tore, liens sphériques de tous les sommets, degrés 12/14/15/16, faces
pentagonales/hexagonales et hexagones isolés. Résultat : **121 PASS, 0 FAIL**.
Le fichier compact `certificates/catalogue_topology.npz` permet de refaire
ces vérifications sans télécharger l'archive externe. Cela certifie
l'existence des pavages, pas leurs minima d'énergie.

**Les six bornes supérieures sont confirmées par un second calcul.**
`scripts/audit_certificates.py` réintègre les tétraèdres avec des formules
barycentriques rationnelles, contrôle la topologie et borne le transport
affine vers les boîtes exactes. Les six valeurs à neuf décimales sont
inchangées, y compris pour μ dans une boîte hexagonale contenant √3.
Les cinq écarts aux valeurs arXiv dépassent aussi l'arrondi de publication
de 0,5 × 10^-6. Les tests de rejet de maillages altérés et d'invariance
géométrique passent (6 tests).

Les articles 2010 et PNAS fournis ensuite par l'utilisateur sont archivés
localement. Le PDF PNAS reçu compte six pages ; ses tableaux détaillés sont
dans le SI Appendix, encore distinct du fichier reçu.

**Exploration sans symétrie imposée.** Deux perturbations indépendantes des
sommets (1,5 % de R₀) et de la boîte (environ 1 %) pour σ, E11, H et μ ont été
relaxées jusqu'au niveau 2. Les 24 étapes restent géométriquement valides au
test flottant. Les deux valeurs finales d'une même structure diffèrent
d'au plus 1,6 × 10^-9 ; σ reste la plus basse de ces essais. Résultats dans
`data/symmetry_breaking_restarts.csv`. Ce test n'est pas une preuve globale.

**Travail analytique.** La combinaison de l'inégalité isopérimétrique, du
minimum du moment quadratique à volume fixé, et de Hölder donne F ≥ 1 pour
toute partition admissible. La borne est commune à toutes les topologies et
ne suffit donc pas à séparer σ de ses concurrentes. Les bornes obtenues à ce
stade pour l'infimum sur le catalogue sont [1, 1,068575909]. La preuve et
l'obstacle à un renforcement purement cellule par cellule sont en rédaction.

## 2026-10-24 — 4. Preuves, dégénérescences et mise au clair du résultat

Le manuscrit `paper.tex` donne maintenant les preuves de la borne universelle,
du critère d'immersion par le degré et du transport affine des certificats.
Il précise les hypothèses de variété et traite aussi le squelette des
triangulations, que le seul contrôle de volume ne suffit pas à justifier.

Un calcul analytique supplémentaire porte sur la famille A15 à faces planes
de l'équation (A16) du préprint de Reddy. Des coefficients de Bernstein
strictement positifs, bornés par Arb à 192 bits sur six intervalles couvrant
[0,2], puis le signe explicite de la dérivée sur [2,+∞[, prouvent
**F(c) > 1,0699 pour tout c ≥ 0**. Le témoin σ certifié est donc meilleur que
tous les membres admissibles de cette famille particulière. Cette preuve
ne couvre pas les partitions A15 courbes. Script : `check_a15_family.py` ;
certificat et journal d'exécution sauvegardés.

Quand c tend vers zéro, deux cellules Z12 de cette famille disparaissent
alors que F tend vers **1,157131807984…**. Une énergie bornée n'exclut donc
pas à elle seule la disparition de cellules. Cette limite précise est
au-dessus de σ ; elle ne résout pas les autres dégénérescences. Le manuscrit
énonce aussi les conditions stationnaires et une forme suffisante de
certificat inférieur global par champs de divergence ; les champs nécessaires
près du niveau σ n'ont pas été construits.

Les README, `result.md`, `references.md` et `NEEDS.md` sont harmonisés avec
les sources effectivement reçues. Le nom de colonne CSV
`F_maxlevel_upper_bound` devient `F_maxlevel_mesh_energy` sans modifier les
valeurs numériques. `verification.md` est explicitement un audit par le
même agent, pas une revue indépendante. La résolution globale demandée
reste non obtenue : aucune borne inférieure ne sépare encore les infima
des topologies concurrentes de la borne supérieure de σ.

## 2026-10-24 — 5. Supplément PNAS retrouvé et contrôle des valeurs publiées

Le SI Appendix a finalement été téléchargé par l'API publique Europe PMC
(`PMC6187128/supplementaryFiles`) ; aucun fichier supplémentaire n'est
nécessaire de la part de l'utilisateur pour cette comparaison. Le PDF compte
17 pages et son SHA-256 est enregistré dans `data/pnas_si_crosscheck.json`.

Lecture et contrôle visuel des pages 4–5 (formule A15, méthode) et 6–12
(figures S3–S14). Les douze énergies DFM non contraintes et les rapports de
boîte correspondent aux entrées du préprint utilisées par le projet.
**Les cinq améliorations certifiées restent donc valables vis-à-vis du SI
publié.** Le script `check_source_comparison.py` vérifie la cohérence de la
transcription manuelle avec les données historiques et les seuils des
certificats. Il ne remplace pas la lecture visuelle de la source.

Les légendes signalent que les calculs hexagonaux utilisent dans Surface
Evolver une période orthorhombique doublée. Un témoin hexagonal se répète
dans cette période sans changer F : A, J et V doublent tous dans l'identité
F³ = 20 A² J / (27 V³). Ce détail est maintenant explicite dans le manuscrit.

L'équation S16 conserve la formule A15 de l'équation A16 du préprint.
Les documents de synthèse et le suivi des sources sont mis à jour ; les
mentions antérieures de SI manquant ci-dessus décrivent les étapes passées.

## 2026-10-24 — 6. État des livrables et contrôle final

Le manuscrit de **neuf pages** est compilé avec le préambule du dépôt, sans
débordement signalé par LaTeX. Les neuf pages rendues ont été inspectées ;
une coupure mal placée d'énoncé a été corrigée. Le contrôle structurel par
pikepdf/qpdf ne signale aucune erreur ni avertissement. Les tests numériques,
certificats et contrôles de sources sont sauvegardés avec leurs journaux.
`git diff --check` ne signale pas d'erreur.

Livrables principaux : `paper.pdf`, `paper.tex`, `verification.md`, les
certificats autonomes du catalogue, les nouveaux vérificateurs et ce journal.
Les PDF tiers restent dans le répertoire local ignoré par Git.

**Résultat à ce stade : partiel.** Les 121 réalisations, les six témoins
d'énergie et la borne restreinte A15 sont certifiés par les calculs décrits.
Le classement global des infima n'est ni prouvé ni réfuté. Il manque une
borne inférieure assez forte sur les classes concurrentes, incluant leurs
limites dégénérées ; aucun article supplémentaire identifié ne ferme ce
point. Aucun sous-agent n'a été utilisé.
