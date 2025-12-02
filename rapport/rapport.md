*Joachim Larrouy (M1 ARIAS Recherche), Nathan Rissot (M1 ARIAS GPEx)*

## Roy-Floyd-Warshall

L'objectif de l'algorithme de Roy-Floyd-Warshall (parfois simplement appelé Floyd dans la suite de ce rapport) est de construire en temps $\Theta (|V|^3)$ `mat_distance` la matrice des plus courtes distances définie comme

$$
\forall i,j \in V, \mathtt{mat\_distance}[i][j] = \min(\omega (i,*,j))
$$
(avec $\omega (i,*,j)$ le plus court chemin de $i$ à $j$ dans le graphe) Cet algorithme prend en entrée la matrice d'adjacence `mat_adjacence` préparée du graphe. On définit `mat_adjacence` comme

$$
\forall i,j \in V, \mathtt{mat\_adjacence}[i][j] = 
\begin{cases} 
    0 \text{ si } i = j \\
    \omega(i,j) \text{ si } (i,j) \in E \\
    \infin \text{ sinon} 
\end{cases} 
$$

La version séquentielle de cet algorithme consiste en une triple boucle $l$,$i$,$j$ de $0$ à $|V|$ ou, pour chaque combinaison de $l$,$i$,$j$ :

$$
\mathtt{mat\_distance}[i][j] = \min (\mathtt{mat\_distance}[i][j], \mathtt{mat\_distance}[i][l] + \mathtt{mat\_distance}[l][j])
$$

La parallélisation de cet algo repose sur le découpage en $x$ blocs de taille $b \times b$ (avec $b = \frac{ |V| }{\sqrt{x}}$).

En effet, à chaque itération de la boucle, un bloc n'a besoin que d'une ligne de données reçue de sa colonne de blocs, et d'une colonne de donnée reçue de sa ligne de blocs. En partageant donc les blocs entre $x$ coeurs, et en les faisant partager les lignes et colonnes de donnée utiles, on peut effectuer les opérations de comparaisons en simultané, ce qui accélere l'exécution.

## PAM

### Presentation général.

L'algorithme PAM est un algorithme ayant pour but de déterminer dans notre ensemble de noeud (noté $V$) un sous ensemble $K$ (que l'on va nommer ensemble de medoïdes) telle que le cout de notre sous ensemble $K$ soit le plus petit possible. Le cout d'un sous ensemble est définit de la manière suivante :

$$
cout = \sum_{n = 0}^{nb\_nodes} \min (\omega(node_n,k_{i}))
$$

- avec $k_i$ un noeud appartenant au sous ensemble $K$.
- $node_n$ le $n$-ième noeud.
- $\omega(n_1,n_2)$ la distance entre deux noeuds au seins de notre matrice des distances.

Pour ce faire, l'algorithme choisit effectue les opération suivantes.

- initialiser notre ensemble de medoïde avec un tirage de noeud aléatoire.
- Pour chaque noeud $n$ appartenent à $( \{V \\ K\} )$ On calcule le cout du remplacement de chaque medoïde par $n$. Si le cout le plus faible entre toutes nos substitution $(swap(n,k_i) = Cout_{min} )$ est plus faible que le cout de l'ensemble $K$ à l'iteration courante, alors on substitue $k_i$ par $n$.

Une fois tout les noeuds evalué, on considère avoir un sous ensemble $K$ représentant une solution satisfaisante à notre recherche de cout minimum.

### Notre heuristique
Le cout de l'algorithme glouton effectuant le calcul des meilleurs médoïdes jusqu'a stabilisation ayant un cout $O(n^2*k^2)$ nous avons choisit d'appliquer une initialisation avec l'heuristique suivant.

Lors de l'initialisation, plutot que de partir de noeud aléatoire, une approche pourrait être de découper en ligne notre matrice des distances, et de déterminer pour chaque fragment un ensemble de noeud qui minimise le cout local a ce fragment.

Puis on choisit comme noeud de départ pour les médoïdes, les noeuds ayant était choisis le plus de fois parmis comme faisant partis des noeuds minimisant les couts locaux.

Comme le calcul du cout favorise les noeuds présent au sein du fragment et que chaque fragment à une probabilité $\frac{1}{nb\_fragment}$ de posséder un médoïde $k_i$. Il convient de retirer du choix possible les noeuds présent dans notre fragment.

![fig x. Représentation du graphe global](./asset/image_ensemble.png)

![fig x. Représentation du fragment local](./asset/image_local.png)


Le cout étant de l'ordre de $n^2$, fragmenter suffisament la matrice permet de faire un précalcul relativement rapide.

> **Remarque** : Dans le cas ou la taille des données reste trop grande pour effectuer ce prétraitement, il est aussi possible de n'effectuer celui-ci que sur un échaintillon de notre donnée.

## Application à des données de séquences d’ARN

Pour l'applications des algorithmes à ce nouveau type de donnée, la seule adaptation que nous avons du faire à été de programmer un nouveau parseur pour extraire les données des fichiers sources, et de mettre en place la génération de la matrice d'adjacence sur ces données.

### Parsing du fichier

Les fichiers aillant une structure assez simple, le parseur en lui même (contenu dans `readArnFromFile`) consiste simplement à placer et à aligner les chaines de caractères extraite du fichier dans une matrice de taille $nb\_nodes \times ||séquence\ arn||$. Cette opération n'est pas parallélisable, toutes les données etant stockés dans le même fichier.

### Construction de la matrice d'adjacence

La construction de la matrice d'adjacence en revanche est parallélisable:

1. ROOT partages les séquences d'ARN avec un `Scatterv` de sorte que chaque coeur dispose de $(\frac{nb\_nodes}{nprocs}) (\pm 1)$ séquences.
2. Les coeurs `Bcast` leurs fragment à tout les autres coeurs.
3. Les coeur utilisent les données reçue pour calculer la distance de hamming entre chaque séquences reçue et chaque séquence dans leurs propre fragment, puis notent le résultat dans leurs matrice.
4. ROOT `Reduce` avec une somme les matrices (initialisées à 0) de chaque coeur.

> **Remarque** :
> Ce coeur est techniquement inefficace car la matrice est symétrique, ce qui signifie que chaque séquence est comparé avec chaque autre séquence deux fois. Pour regler cette inneficacité il faudrait pouvoir découper la charge de travail de manière triangulaire.

> **Remarque 2** : 
> Dans un soucis de simplicité d'écriture, chaque coeur réserve en mémoire une matrice complète de taille $nb\_nodes \times nb\_nodes$. Cette facilité à été laissée comme amélioration potentielle au moment de l'écriture, et est toujours dans le code aujourd'hui par manque de temps.

Une fois la matrice d'adjacence obtenue, on peut procéder au découpage par blocs de celles-ci pour l'algorithme de Roy-Floyd-Warshall afin d'obtenir la matrice des distances des plus courts chemins comme précedemment, puis au découpage par ligne afin d'utiliser PAM pour obtenir les $k$-médoïdes.