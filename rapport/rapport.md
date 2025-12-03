# Programmation Parallèle - Projet 1 - MPI
*Joachim Larrouy (M1 ARIAS Recherche) - Nathan Rissot (M1 ARIAS GPEx)*

L'objectif de ce projet est, dans un premier temps, de construire une version parallèles des algorithmes de Roy-Floyd-Warshall, puis PAM, afin de pouvoir identifier $k$ communautées centrées autour de $k$-médoïdes.

Puis, dans un second temps, d'adapter le processus de traitements à des données sous la forme de séquences d'ARN

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

![*fig 1. Exemple de partage de la matrice par bloc*](./asset/bloc_share_floyd_schéma_tableau.jpg)

## PAM

### Presentation général.

L'algorithme PAM est un algorithme ayant pour but de déterminer dans notre ensemble de noeud (noté $V$) un sous ensemble $K$ (que l'on va nommer ensemble de medoïdes) telle que le cout de notre sous ensemble $K$ soit le plus petit possible. Le cout d'un sous ensemble est définit de la manière suivante :

$$
\text{cout} = \sum_{n = 0}^{nb\_nodes} \min (\omega(node_n,k_{i}))
$$

- avec $k_i$ un noeud appartenant au sous ensemble $K$.
- $node_n$ le $n$-ième noeud.
- $\omega(n_1,n_2)$ la distance entre deux noeuds au seins de notre matrice des distances.

Pour ce faire, l'algorithme choisit effectue les opération suivantes.

- initialiser notre ensemble de medoïde avec un tirage de noeud aléatoire.
- Pour chaque noeud $n$ appartenent à $( \{V \\ K\} )$ On calcule le cout du remplacement de chaque medoïde par $n$. Si le cout le plus faible entre toutes nos substitution $(\text{swap}(n,k_i) = \text{cout}_{\min} )$ est plus faible que le cout de l'ensemble $K$ à l'iteration courante, alors on substitue $k_i$ par $n$.

Une fois tout les noeuds evalué, on considère avoir un sous ensemble $K$ représentant une solution satisfaisante à notre recherche de cout minimum.

### Notre heuristique
Le cout de l'algorithme glouton effectuant le calcul des meilleurs médoïdes jusqu'a stabilisation ayant un cout $O(n^2*k^2)$ nous avons choisit d'appliquer une initialisation avec l'heuristique suivant.

Lors de l'initialisation, plutot que de partir de noeud aléatoire, une approche pourrait être de découper en ligne notre matrice des distances, et de déterminer pour chaque fragment un ensemble de noeud qui minimise le cout local a ce fragment.

Puis on choisit comme noeud de départ pour les médoïdes, les noeuds ayant était choisis le plus de fois parmis comme faisant partis des noeuds minimisant les couts locaux.

Comme le calcul du cout favorise les noeuds présent au sein du fragment et que chaque fragment à une probabilité $\frac{1}{nb\_fragment}$ de posséder un médoïde $k_i$. Il convient de retirer du choix possible les noeuds présent dans notre fragment.

<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1em;"><figure><img src="./asset/image_ensemble.png" alt="fig 2. Représentation du graphe "/><figcaption aria-hidden="true"><i>fig 2. Représentation du graphe global</i></figcaption></figure><figure><img src="./asset/image_local.png" alt="fig 3. Représentation du fragment "/><figcaption aria-hidden="true"><i>fig 3. Représentation du fragment local</i></figcaption></figure></div>

Le cout étant de l'ordre de $n^2$, fragmenter suffisament la matrice permet de faire un précalcul relativement rapide.

> **Remarque** : 
> Dans le cas ou la taille des données reste trop grande pour effectuer ce prétraitement, il est aussi possible de n'effectuer celui-ci que sur un échaintillon de notre donnée.

## Application à des données de séquences d’ARN

Pour l'applications des algorithmes à ce nouveau type de donnée, la seule adaptation que nous avons du faire à été de programmer un nouveau parseur pour extraire les données des fichiers sources, et de mettre en place la génération de la matrice d'adjacence sur ces données.

### Parsing du fichier

Les fichiers aillant une structure assez simple, le parseur en lui même (contenu dans `readArnFromFile`) consiste simplement à placer et à aligner les chaines de caractères extraite du fichier dans une matrice de taille $nb\_nodes \times ||séquence\ arn||$. Cette opération n'est pas parallélisable, toutes les données etant stockés dans le même fichier.

### Construction de la matrice d'adjacence

La construction de la matrice d'adjacence en revanche est parallélisable:

1. ROOT partages les séquences d'ARN avec un `Scatterv` de sorte que chaque coeur dispose de $\left(\frac{nb\_nodes}{nprocs}\right) (\pm 1)$ séquences.
2. Les coeurs `Bcast` leurs fragment à tout les autres coeurs.
3. Les coeur utilisent les données reçue pour calculer la distance de hamming entre chaque séquences reçue et chaque séquence dans leurs propre fragment, puis notent le résultat dans leurs matrice.
4. ROOT `Reduce` avec une somme les matrices (initialisées à 0) de chaque coeur.

> **Remarque** :
> Ce coeur est techniquement inefficace car la matrice est symétrique, ce qui signifie que chaque séquence est comparé avec chaque autre séquence deux fois. Pour regler cette inneficacité il faudrait pouvoir découper la charge de travail de manière triangulaire.

> **Remarque** : 
> Dans un soucis de simplicité d'écriture, chaque coeur réserve en mémoire une matrice complète de taille $nb\_nodes \times nb\_nodes$. Cette facilité à été laissée comme amélioration potentielle au moment de l'écriture, et est toujours dans le code aujourd'hui par manque de temps.

Une fois la matrice d'adjacence obtenue, on peut procéder au découpage par blocs de celles-ci pour l'algorithme de Roy-Floyd-Warshall afin d'obtenir la matrice des distances des plus courts chemins comme précedemment, puis au découpage par ligne afin d'utiliser PAM pour obtenir les $k$-médoïdes.

## Benchmarks & Analyses

> **Remarque** :
> Nous n'avons effectué des Bench que sur les séquences d'ARN car les temps de calculs sur les exemples de graphes donnés etait trop court, et les différences de temps venaient plus de facteurs extérieurs (temps d'initialisation de MPI, opération OS) que du code lui même.

Les bench ont été effectués sur des TODO:TYPEOFMACHINE diposant de 6 coeurs physiques. Les test sur 16 processeurs ont donc été lancé avec l'option `--oversubscribe`. 

De plus, nous n'avons pas lancé de test avec 9 coeurs, car les tailles des fichiers de séquences d'ARN (500 et 2k) ne sont pas divisible par $\sqrt{9} (=3)$ et donc ne respectent pas les hypothèses favorables dont dépend notre implémentation de l'algorithmes de Roy-Floyd-Warshall

![*fig 4. Performances moyennes (en ms) des étapes de génération de la matrice, de l'executions de l'algorithmes de Roy-Floyd-Warshall, de PAM et du total pour 500 noeuds selon le nombre de coeurs*](asset/500_nodes.svg)

![*fig 5. Performances moyennes (en ms) des étapes de génération de la matrice, de l'executions de l'algorithmes de Roy-Floyd-Warshall, de PAM et du total pour 2K noeuds selon le nombre de coeurs*](asset/2K_nodes.svg)

Sur ces graphes des temps d'execution moyens, on remarque que pour presque toutes les étapes, le temps de calcul diminue au fur et à mesure que l'on augmente le nomnbre de processus (CF, fig 8. les courbes d'accélérations). 

Cepandant on remarque que pour PAM, le temps de calcul augmente légerement entre la version séquentielle et la version parallélisée. Cela est dù au fait que l'implémentation du premier calcul qui nous permet d'avoir une initialisation fiable des premiers candidats ne fonctionne pas avec un seul processeur, et à donc été remplacé par un simple choix aléatoire.

![*fig 6. Courbe des temps d'execution cumulés, et pourcentage du temps total pour 500 noeuds selon le nombre de coeurs*](asset/temps_cumulés_500.svg)

Comme on peut le voir sur le graphe des courbes cumulés, ce changement implique qu'une plus grande part du temps d'execution total est dédié a cette étape, mais on peut aussi constater que les résultats obtenus par les executions de l'algorithme avec le pré-traitement heuristique sont bien meilleurs

![*fig 7. Comparaison des temps d'exécutions, et de la quelité des résultats de PAM avec et sans l'heuristique*](asset/TODO.svg)

![*fig 8. Courbe d'accélération de l'execution du programme*](asset/acceleration_2k.svg)

> **Calcul de l'accelération** : 
> $$\text{speedup} = \frac{T_{seq}}{T_{nprocs}} \times 100$$