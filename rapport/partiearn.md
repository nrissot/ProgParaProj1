## Application à des données de séquences d’ARN

Pour l'applications des algorithmes à ce nouveau type de donnée, la seule adaptation que nous avons du faire à été de programmer un nouveau parseur pour extraire les données des fichiers sources, et de mettre en place la génération de la matrice d'adjacence sur ces données.

### Parsing du fichier

Les fichiers aillant une structure assez simple, le parseur en lui même (contenu dans `readArnFromFile`) consiste simplement à placer et à aligner les chaines de caractères extraite du fichier dans une matrice de taille $nb\_nodes \times ||séquence\ arn||$. Cette opération n'est pas parallélisable, toutes les données etant stockés dans le même fichier.

### Construction de la matrice d'adjacence

La construction de la matrice d'adjacence en revanche est parallélisable:

1. ROOT partages les séquences d'ARN avec un `Scatterv` de sorte que chaque processus dispose de $(\frac{nb\_nodes}{nprocs}) (\pm 1)$ séquences.
2. Les processus `Bcast` leurs fragment à tout les autres processus.
3. Les processus utilisent les données reçue pour calculer la distance de hamming entre chaque séquences reçue et chaque séquence dans leurs propre fragment, puis notent le résultat dans leurs matrice.
4. ROOT `Reduce` avec une somme les matrices (initialisées à 0) de chaque processus.

> **Remarque** :
> Ce processus est techniquement inefficace car la matrice est symétrique, ce qui signifie que chaque séquence est comparé avec chaque autre séquence deux fois. Pour regler cette inneficacité il faudrait pouvoir découper la charge de travail de manière triangulaire.

> **Remarque 2** : 
> Dans un soucis de simplicité d'écriture, chaque processus réserve en mémoire une matrice complète de taille $nb\_nodes \times nb\_nodes$. Cette facilité à été laissée comme amélioration potentielle au moment de l'écriture, et est toujours dans le code aujourd'hui par manque de temps.

Une fois la matrice d'adjacence obtenue, on peut procéder au découpage par blocs de celles-ci pour l'algorithme de Roy-Floyd-Warshall afin d'obtenir la matrice des distances des plus courts chemins comme précedemment, puis au découpage par ligne afin d'utiliser PAM pour obtenir les $k$-médoïdes.