## PAM

### Presentation général.

L'algorithme PAM est un algorithme ayant pour but de déterminer dans notre ensemble de noeud (noté $V$) un sous ensemble $K$ (que l'on va nommer ensemble de medoïdes) telle que le cout de notre sous ensemble $K$ soit le plus petit possible. Le cout d'un sous ensemble est définit de la manière suivante :

$$cout = \sum_{n = 0}^{nb\_nodes} min(distance(node_n,k_{i}))$$

- avec $k_i$ un noeud appartenant au sous ensemble $K$.
- $node_n$ le $n$-ième noeud.
- $distance(n_1,n_2)$ la distance entre deux noeuds au seins de notre matrice des distances.

Pour ce faire, l'algorithme choisit effectue les opération suivantes.
- initialiser notre ensemble de medoïde avec un tirage de noeud aléatoire.
- Pour chaque noeud $n$ appartenent à $( \{V \\ K\} )$ On calcule le cout du remplacement de chaque medoïde par $n$. Si le cout le plus faible entre toutes nos substitution $(swap(n,k_i) = Cout_{min} )$ est plus faible que le cout de l'ensemble $K$ à l'iteration courrante, alors on substitue $k_i$ par $n$.

Une fois tout les noeuds evalué, on considère avoir un sous ensemble $K$ représente une solution satisfaisante à notre recherche de cout minimum.

### Notre heuristique
Le cout de l'algorithme glouton effectuant le calcul des meilleurs médoïdes jusqu'a stabilisation ayant un cout $O(n^2*k^2)$ nous avons choisit d'appliquer une initialisation avec l'heuristique suivant.

Lors de l'initialisation, plutot que de partir de noeud aléatoire, une approche pourrait être de découper en ligne notre matrice des distances, et de déterminer pour chaque fragment un ensemble de noeud qui minimise le cout local a ce fragment.

Puis on choisit comme noeud de départ pour les médoïdes, les noeuds ayant était choisis le plus de fois parmis comme faisant partis des noeuds minimisant les couts locaux.

Comme le calcul du cout favorise les noeuds présent au sein du fragment et que chaque fragment à une probabilité $\frac{1}{nb\_fragment}$ de posséder un médoïde $k_i$. Il convient de retirer du choix possible les noeuds présent dans notre fragment.

![fig x. Représentation du graphe global](./asset/image_ensemble.png)

![fig x. Représentation du fragment local](./asset/image_local.png)


Le cout étant de l'ordre de $n^2$, fragmenter suffisament la matrice permet de faire un précalcul relativement rapide.

Remarque : Dans le cas ou la taille des données reste trop grande pour effectuer ce prétraitement, il est aussi possible de n'effectuer celui-ci que sur un échaintillon de notre donnée.