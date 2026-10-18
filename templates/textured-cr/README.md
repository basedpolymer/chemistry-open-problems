# Gabarit TeXtured des comptes rendus

Compilation : `python3 .claude/scripts/cr_pdf.py "<CR.md>"`. Le PDF voisin du Markdown est l'aperçu des CR ;
Markdown Preview Enhanced reste réservé aux notes. `--suivi` recompile après modification du Markdown, des figures
ou du gabarit. `--garder` conserve les fichiers LaTeX et le journal ; `--sous-titre` remplace le sous-titre du document.

## Comptes rendus structurés

Le mode structuré est activé dans l'en-tête YAML du Markdown ; les documents sans cette métadonnée gardent
leur mise en page historique. Le sous-titre du YAML est utilisé par défaut.

```yaml
---
cr-structured: true
subtitle: "Groupe F · Séance … · Encadrante …"
---
```

Un seul titre `#`. Les parties `##` et les sous-parties `###` n'ont pas de numéros saisis : LaTeX les calcule.
Ajouter des identifiants stables aux titres ; les parties liminaires et finales portent `.unnumbered`.

```markdown
## Introduction {#sec:introduction .unnumbered}
## Comportement mécanique en flexion {#sec:flexion}
### Méthode {#sec:methode-flexion}
```

Renvois : `\CRsection{sec:methode-flexion}` affiche le numéro de section cliquable ;
`\CRref{tableau}{tab:flexion}` et `\CRref{figure}{fig:tp3F_flexion}` affichent le nom et le numéro.
Le filtre `cr-structure.lua` ajoute aux figures une cible `fig:<nom du fichier sans extension>`, sans modifier
leur légende, leur taille, leur fichier ou leurs réglages de placement. Un identifiant explicite existant est conservé.
Les renvois aux pages et figures d'un support externe restent du texte, pour éviter de les confondre avec ceux du CR.

## Tableaux, nombres et unités

`CRtable` garde ensemble la légende et le tableau court. `tabularray` fournit les colonnes ; `S` aligne les nombres
avec `siunitx`. Les en-têtes utilisent `guard` pour ne pas être interprétés comme des nombres. Les grandes tables
multipages continuent à employer `longtable` ou `longtblr`.

```latex
\begin{CRtable}{Modules de flexion (\si{\giga\pascal}).}{tab:flexion}
\begin{tblr}{width=\linewidth, colspec={X[l] S[table-format=1.2(2)]}, row{1}={guard}}
\toprule
Matériau & Module \\
\midrule
PMMA & 2.3 +- 0.1 \\
PMMA choc & 1.70 +- 0.08 \\
\bottomrule
\end{tblr}
\end{CRtable}
```

Dans les commandes `siunitx`, saisir les décimales avec un point ; le PDF affiche une virgule et conserve les
décimales significatives : `\num{1.70 +- 0.08}`, `\SI{1.70 +- 0.08}{\giga\pascal}`,
`\SIrange{7}{11}{\percent}`. `\si{\kilo\joule\per\metre\squared}` affiche l'unité seule.
`\SI` est utilisé pour éviter le conflit de `\qty` avec le paquet `physics` du préambule original.

Les figures, tableaux et équations ont des compteurs distincts. `\CRappendix` avant la première annexe remet
le compteur des tableaux à zéro et produit A1, A2, etc. ; un nouvel appel passe à B1, B2, etc.

## Encadrés et vérification

Les encadrés `CRcalcul` et `CRremarque` ont un titre, un filet gris et aucun numéro. Leur contenu est en LaTeX ;
reprendre un raisonnement existant, sans y introduire de résultat supplémentaire.

```latex
\begin{CRremarque}{Conditions de validité}
Texte rédigé, avec formules en \( \) ou \[ \].
\end{CRremarque}
```

Une annotation de travail qui doit rester hors du PDF peut être conservée dans un bloc brut LaTeX :

````markdown
```{=latex}
% TODO : point de rédaction à reprendre
```
````

Contrôler le PDF, ses liens, les avertissements de compilation, les unités et les chiffres significatifs.
Les environnements natifs `Todo` de TeXtured restent visibles même avec `\WIPfalse`.
La bibliographie reste la liste écrite dans le Markdown ; aucun traitement bibliographique n'est ajouté.
