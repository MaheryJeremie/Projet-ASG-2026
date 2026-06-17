# Pipeline d'Assemblage de Génome De Novo

## 1. Comment tester l'application ?

Pour lancer et tester l'application sur votre machine :

1. Ouvrez un terminal dans le dossier du projet :

```bash
...\Projet ASG-2026
```

2. Lancez la commande suivante :

```bash
streamlit run Home.py
```

3. Votre navigateur internet s'ouvrira automatiquement à l'adresse :

```
http://localhost:8501
```

sur la page d'accueil de l'application (`Home.py`).

4. Utilisez le menu latéral gauche pour naviguer et tester les trois modules (Lots 1, 2 et 3).

5. Dans chaque module, vous pouvez :
   - Cliquer sur **"Générer le Toy Dataset"**
   - Utiliser les exemples intégrés

afin d'observer les algorithmes s'exécuter en temps réel avec graphiques et statistiques.

---

## 2. À quoi sert cette application ?

Cette application implémente un **pipeline d'assemblage de génome *De Novo***.

### Le problème biologique

Les machines de séquençage ADN modernes ne peuvent pas lire un génome complet (contenant des millions voire des milliards de nucléotides A, C, G et T) en une seule opération.

Elles produisent donc des millions de petits fragments aléatoires appelés **reads**, généralement longs de 50 à 150 nucléotides.

### Objectif

L'objectif est de reconstituer la séquence ADN originale à partir de ces fragments désordonnés, comme un immense puzzle.

La séquence reconstruite est appelée un **contig**.

Cette approche est dite **De Novo** car elle ne s'appuie sur aucun génome de référence.

### Applications

- Identification de nouveaux virus
- Étude d'organismes inconnus
- Recherche en génomique
- Analyse de biodiversité

---

# 3. Fonctionnement détaillé des trois lots

## Lot 1 : Ingestion & Contrôle Qualité

**Fichiers concernés :**

- `1_Ingestion_qualite.py`
- `kmer.py`

### Objectif

Lire les fichiers de séquençage (FASTQ) et éliminer les lectures de mauvaise qualité.

### Étapes

#### Lecture des reads

Les reads sont extraits depuis les fichiers FASTQ avec leurs scores de qualité associés.

#### Filtrage qualité

Les reads présentant des scores de qualité Phred trop faibles sont supprimés.

#### Analyse par k-mers

Chaque read est découpé en sous-séquences de taille fixe `k`, appelées **k-mers**.

Exemple :

```
Read : ACGTACG
k = 3

ACG
CGT
GTA
TAC
ACG
```

#### Histogramme des fréquences

L'application compte le nombre d'occurrences de chaque k-mer.

On distingue :

- **Hapax** : k-mers observés une seule fois (souvent dus à des erreurs de lecture)
- **K-mers solides** : k-mers apparaissant fréquemment et représentant le véritable signal biologique

L'application calcule automatiquement un seuil permettant de séparer les erreurs du signal réel.

---

## Lot 2 : Alignement par Programmation Dynamique (LCS)

**Fichiers concernés :**

- `2_Alignement_lcs.py`
- `alignment.py`

### Objectif

Déterminer si deux reads proviennent d'une même région du génome en recherchant un chevauchement.

### Algorithme utilisé : LCS

LCS signifie :

**Longest Common Subsequence**  
(**Plus Longue Sous-Séquence Commune**)

L'algorithme repose sur la programmation dynamique afin de trouver les similarités entre deux séquences même en présence :

- de mutations
- d'insertions
- de suppressions (délétions)

### Visualisations proposées

L'application affiche :

- La matrice de programmation dynamique
- L'alignement reconstruit

Exemple :

```text
ACGT-ACG
|||  |||
ACGTTACG
```

Les symboles :

- `|` : caractères identiques
- `-` : insertion ou délétion

---

## Lot 3 : Assemblage Memory-Efficient

### Bloom Filter & Graphe de de Bruijn

**Fichiers concernés :**

- `3_Assemblage_bloom.py`
- `assembly.py`

### Objectif

Assembler les reads en longues séquences (contigs) tout en limitant fortement l'utilisation de la mémoire.

---

### Filtre de Bloom

Stocker plusieurs millions de k-mers dans une structure classique (dictionnaire ou ensemble) consomme énormément de RAM.

Pour résoudre ce problème, l'application utilise un **Filtre de Bloom**.

#### Principe

Le filtre de Bloom :

- utilise un tableau de bits
- applique plusieurs fonctions de hachage (`mmh3`)
- stocke uniquement l'information de présence probable

### Avantages

- Très faible consommation mémoire
- Facteur de réduction de 10 à 20 fois par rapport à un stockage classique
- Vérification rapide de la présence d'un k-mer

---

### Graphe de de Bruijn implicite

Traditionnellement, l'assemblage repose sur la construction d'un immense graphe reliant les k-mers qui se chevauchent.

Dans cette application, le graphe n'est jamais stocké explicitement.

### Méthode utilisée

À partir d'un k-mer de départ (*seed*) :

1. On considère les quatre extensions possibles :

```
A
C
G
T
```

2. Pour chaque extension candidate, on interroge le filtre de Bloom :

> « Ce nouveau k-mer existe-t-il ? »

3. Si la réponse est positive :
   - l'extension est ajoutée à la séquence
   - l'assemblage continue

4. L'algorithme s'arrête lorsqu'il rencontre :
   - une impasse
   - un embranchement (bifurcation)

---

### Résultat

Le processus produit des **contigs**, c'est-à-dire de longues séquences ADN reconstituées à partir des reads initiaux.

Ces contigs représentent les fragments les plus longs du génome reconstruits par le pipeline.

---

# Résumé du Pipeline

```text
FASTQ
  │
  ▼
Lot 1 : Contrôle Qualité
  │
  ▼
Extraction des k-mers
  │
  ▼
Filtrage des k-mers solides
  │
  ▼
Lot 2 : Alignement LCS
  │
  ▼
Détection des chevauchements
  │
  ▼
Lot 3 : Bloom Filter
  │
  ▼
Graphe de de Bruijn implicite
  │
  ▼
Assemblage
  │
  ▼
Contigs (génome reconstruit)
```
