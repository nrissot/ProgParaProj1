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

La parallélisation de cet algo repose sur le découpage en $x$ blocs de taille $b \times b$ (avec $b = \frac{ |V| }{\sqrt{x}}$) et l'application de l'algorithme en paralléle sur chacun de ces blocs. En effet, seules les valeurs provenant d'une seule des lignes et d'une seule des colonnes de la matrice des blocs n'est nécessaire pour qu'un bloc calcule les nouveaux plus courts chemins à une itération $k$ donnée.

