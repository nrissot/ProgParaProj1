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

