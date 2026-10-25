# Journal de reprise — Frank–Kasper

## 2026-10-24 — 1. État initial et objectif

Reprise demandée par le propriétaire : viser une résolution complète, sans
sous-agents, et documenter progressivement l'avancement. Travail effectué par
GPT-6-Astra (OpenAI ; effort xhigh). Aucune vérification par un autre agent
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

## 2026-10-25 — 7. Reprise : écart uniforme au modèle sphérique

La demande de résolution complète est reprise au niveau des minorants, sans
sous-agents. Répéter les optimisations locales ne fermerait pas ce verrou.
L'analyse porte sur une éventuelle suite de partitions dont F tendrait vers 1,
en autorisant toutes les boîtes et des volumes positifs arbitrairement petits.

Une normalisation des rayons équivalents r_i par la moyenne harmonique,
pondérée par les volumes, donne une estimation simultanée de leur dispersion
et du volume manquant aux boules centrées aux centroïdes. Les boules des
cellules satisfaisant deux seuils explicites peuvent ensuite être rétrécies
en un empilement de boules égales. La borne élémentaire de Blichfeldt, dont
une preuve complète est incluse, ou le théorème de Kepler, interdisent que
cet empilement occupe presque tout l'espace.

Le script `check_uniform_gap.py` vérifie uniquement les constantes rationnelles
de cet argument. Il ne doit pas être confondu avec une formalisation de la
preuve géométrique ni avec une nouvelle vérification de Flyspeck. Résultat
et contrôle du manuscrit à consigner après exécution. Cette piste vise une
borne universelle très faible, pas encore une séparation entre topologies.

**Contrôle des constantes : PASS.** La preuve élémentaire donne
F > 1,000000004 ; l'utilisation du théorème de Kepler donne
**F > 1,00000025**. Les inégalités entre constantes sont vérifiées avec
`Fraction`, sans arrondi flottant. L'infimum du catalogue appartient donc
à [1,00000025 ; 1,068575909]. Le gain numérique est minuscule ; il ferme
le cas d'une suite tendant vers 1, mais pas les dégénérescences proches de σ.
L'énoncé de Kepler a été lu à la source (DOI 10.1017/fmp.2017.1, section 3
du préprint) ; sa formalisation n'a pas été réexécutée. Le manuscrit contient
la preuve complète de la réduction DFM et de la version élémentaire.

## 2026-10-25 — 8. Résolution globale d'une sous-famille : transformations affines

Pour une partition **fixée**, toutes ses déformations affines se paramètrent
par une matrice symétrique positive G = TᵀT. Les facteurs de déterminant
s'annulent dans F³. Après normalisation du moment, minimiser F revient à
minimiser une fonction convexe de G sous une contrainte linéaire.
Une matrice duale positive fournit donc un véritable minorant global sur
cette famille, et non un simple test stationnaire.

Le nouveau script `check_affine_orbits.py` intègre le tenseur des moments en
rationnels et certifie le minorant par Arb à 192 bits. Les six matrices de
métrique sont enregistrées en rationnels ; une seconde vérification peut les
relire sans relancer l'optimiseur. Résultats :

| Maillage fixé | Infimum sur toutes ses déformations affines |
|---|---|
| BCC | [1,071107807 ; 1,071107811] |
| P | [1,068892909 ; 1,068892956] |
| M | [1,068962667 ; 1,068962682] |
| H | [1,068948492 ; 1,068948499] |
| μ | [1,069099093 ; 1,069099122] |
| σ | [1,068575891 ; 1,068575910] |

Ainsi, aucune déformation affine des cinq maillages concurrents sauvegardés
ne bat le témoin σ. **Cela ne couvre pas leurs autres géométries internes.**
Le calcul historique avec relaxation interne de H est d'ailleurs plus bas
que la borne affine du maillage H à boîte fixée : ce sont bien deux domaines
d'optimisation distincts. Trois tests analytiques supplémentaires passent
(cube, tenseur sous cisaillement et translation, rejet de faux certificats).

La comparaison générale des infima des topologies reste non résolue.
Les nouveaux minorants sont uniformes sur leurs domaines explicités ; aucun
ne permet encore de séparer toutes les géométries concurrentes de σ.

## 2026-10-25 — 9. Vérification et portée finale de cette reprise

La seconde exécution des six certificats affines, utilisant seulement les
métriques rationnelles sauvegardées, reproduit les mêmes intervalles :
6 PASS. Les constantes de la borne uniforme passent également ; les trois
nouveaux tests analytiques sont sauvegardés dans les journaux d'exécution.
Le manuscrit compilé compte treize pages. Le contrôle structurel PDF par
pikepdf/qpdf ne signale ni erreur ni avertissement et LaTeX ne signale aucun
débordement. Les pages rendues sont inspectées, avec les références regroupées
sur une page distincte pour éviter une dernière référence isolée.

Les deux README, le manuscrit, les sources et le rapport d'audit reflètent
les nouveaux résultats. La revue indépendante reste à faire ; aucun
sous-agent n'a été utilisé. **La résolution complète n'a pas été obtenue.**
La recherche a produit un écart uniforme au modèle sphérique et des minorants
globaux pour six familles affines. Elle n'a produit ni minorant séparant σ de
toutes les topologies concurrentes, ni contre-exemple à son optimalité.

## 2026-10-25 — 10. Reprise sur les formes internes et les huit entrées difficiles

À la nouvelle demande de résolution complète, deux pistes sont attaquées.
La première utilise des champs de divergence linéaires par cellule : leurs
différences aux interfaces dépendent seulement des incidences périodiques.
Une optimisation convexe fournit ainsi des candidats minorants valables sur
toutes les formes internes d'une topologie. Une normalisation par le
déterminant permet de traiter aussi toutes les boîtes. Les premiers nombres
restent très inférieurs à σ ; ils sont encore exploratoires à cette étape.

La seconde reprend les huit entrées dont le tableau historique ne conservait
que des configurations initiales. L'optimiseur utilise des cônes positifs,
avec des points intérieurs auxiliaires mobiles, afin de préserver la
réalisation géométrique. Le calcul sur Y04 descend numériquement de 1,6201
à environ 1,2545. Les gradients de la barrière sont contrôlés par différences
centrées ; les bornes décimales de ces nouveaux témoins restent à certifier.

Un point de méthode apparaît : un test de cônes négatif au **centroïde** ne
prouve pas l'invalidité d'une cellule. Les mêmes faces peuvent admettre un
autre point intérieur donnant uniquement des cônes positifs. Les calculs
de Y04 présentent précisément ce comportement. Une vérification entière
avec ces points auxiliaires est nécessaire avant toute conclusion de validité.

## 2026-10-25 — 11. Minorants sur toutes les formes, et plafond de cette méthode

Les champs linéaires par cellule fournissent maintenant 121 certificats
entiers, tous vérifiés avec leurs complexes périodiques. Les contraintes
portent seulement sur les interfaces : |z_i − z_j + W s| ≤ 1. Le théorème
donne, pour toutes les formes internes et toutes les boîtes,
F ≥ (6 det(W)/(π n))^(2/9). Trente-quatre minorants améliorent la borne
uniforme précédente. Exemples : A15 ≥ 1,008666655 ; H ≥ 1,008808885 ;
σ ≥ 1,004288649 ; E11 ≥ 1,004321474. Les propositions numériques sont
arrondies et rétractées avant vérification ; la réussite du solveur n'est
pas une hypothèse du certificat. Une entrée a nécessité le repli sur une
proposition géométrique plus faible après échec du solveur.

Une limite de cette piste est désormais démontrée. Les champs réalisent le
complexe de contact par des arêtes de longueur au plus un. L'identité des
volumes orientés et le volume maximal du tétraèdre régulier imposent
det(W) ≤ t/(6√2). Comme t/n ≤ 124/21 dans le catalogue, le nombre fourni
par cette méthode reste inférieur à 1,065254, donc inférieur au témoin σ
1,068575909. **Ce plafond concerne le minorant construit, pas l'énergie
physique.** Raffiner cette seule optimisation ne pourra donc pas fournir
la séparation recherchée avec le témoin actuel. Il faut exploiter plus
finement le couplage entre aire, moments et incidences.

Les quatre tests des certificats linéaires et les deux tests de gradient
et de frontière de la barrière passent. Les six tests préexistants de
l'auditeur passent après ajout de la comparaison directe au catalogue
entier, sans dépendance à l'archive externe.

## 2026-10-25 — 12. Contre-exemple exact au diagnostic nécessaire par centroïde

Le pré-audit entier de Y04 réussit : 42 cellules, 290 faces, 248 sommets
de mousse, identité exacte du complexe périodique et tous les cônes
auxiliaires strictement positifs. L'énergie du témoin appartient à
[1,254466832434 ; 1,254466832435]. Pourtant, les moments rationnels donnent
**184 cônes négatifs au vrai centroïde**. Ce n'est donc pas seulement une
impression numérique : une partition admissible peut échouer à ce test.
La documentation de `cone_check` est corrigée pour décrire une condition
suffisante de validité, jamais une condition nécessaire. Un troisième test
contrôle que le centroïde calculé pour un cube reste le même lorsqu'on
déplace son sommet auxiliaire intérieur.

## 2026-10-25 — 13. Huit reprises achevées, huit témoins certifiés

Les huit calculs de niveau zéro sont terminés : six poids de barrière par
entrée, soit 48 étapes documentées dans `data/barrier_continuation.csv`.
Les huit maillages arrondis sur la grille 2⁻⁴⁴ passent la comparaison exacte
au complexe catalogue, les tests de variété, la positivité des cônes
auxiliaires, la somme exacte des volumes et l'intégration rationnelle.

| Entrée | Borne supérieure certifiée |
|---|---:|
| N55 | 1,102285105 |
| N59 | 1,101112653 |
| N60 | 1,129410812 |
| II:Y02 | 1,103468432 |
| II:Y03 | 1,103507483 |
| II:Y04 | 1,254466833 |
| II:Y05 | 1,254006234 |
| II:Y06 | 1,255816821 |

Y05 et Y06 fournissent aussi des exemples exacts de cellules admissibles
échouant au test au centroïde : 370 et 281 cônes négatifs, respectivement.
Les huit anciennes lignes du tableau historique restent identifiées comme
telles ; le supplément donne les nouveaux résultats sans leur attribuer
une extrapolation ou une optimalité non démontrée. Les gradients finaux
de certains calculs sont encore élevés. Les nouveaux témoins restent tous
au-dessus de σ et n'excluent pas des infima plus bas dans ces classes.

La piste linéaire a été menée jusqu'à sa limitation démontrée. La comparaison
globale reste ouverte : il manque une inégalité exploitant assez fortement
le couplage entre interfaces, moments et incidences pour séparer les infima,
ou un contre-exemple admissible accompagné d'une comparaison rigoureuse.
**La résolution complète n'est toujours pas obtenue.**

## 2026-10-25 — 14. Vérification du supplément et documentation

La vérification sans solveur reproduit les 121 minorants et certifie le
plafond strict 1,065253040. Les huit témoins passent à nouveau à partir
des seuls fichiers entiers sauvegardés ; leurs énergies concordent avec
les 48 étapes enregistrées. Les quatre tests des champs linéaires, trois
tests de barrière/centroïde et six tests existants de l'auditeur passent.
Les huit annotations de diagnostic du CSV historique sont corrigées sans
modifier ses valeurs numériques ; le générateur utilise le même vocabulaire.

Les deux README, le rapport d'audit et le manuscrit sont actualisés. Le
PDF final contient seize pages : compilation sans débordement ni référence
manquante, contrôle structurel pikepdf/qpdf sans anomalie et inspection des
seize pages rendues. Les sources et les données nécessaires à la nouvelle
vérification sont conservées dans le dépôt ; les propositions flottantes
intermédiaires et les articles restent dans les dossiers ignorés.
