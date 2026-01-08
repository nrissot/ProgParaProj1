# `cost_from_candidate_set`

Itérations de la boucle indépendantes les unes des autres, on peut donc remplacer la boucle `for` par un `parallel for`. Il faut bien penser à mettre l'instruction `atomic` avant l'incrémentation pour proteger l'écriture dans la somme.

```c++
int cost_from_candidate_set(int *mat_distance_fragment, int *candidates, int k, int nb_nodes, int nb_lignes_fragment)
{
    int sum = 0;
#pragma omp for
    for (int current_node = 0; current_node < nb_lignes_fragment; current_node++)
    {
        int min_current = INF;
        for (int i = 0; i < k; i++)
        {
            min_current = min(mat_distance_fragment[current_node * nb_nodes + candidates[i]], min_current);
        }
#pragma omp atomic
        sum += min_current;
    }
    return sum;
}
```

# `calculate_cost_fragment`
Très similaire à `cost_from_candidate_set`

```c++
int calculate_cost_fragment(int *medoids, int K, int *mat_distance_fragment, int nb_nodes, int nb_lignes)
{
    int sum = 0;
#pragma omp parallel for 
    for (int i = 0; i < nb_lignes; ++i)
    {
        int min = mat_distance_fragment[i * nb_nodes + medoids[0]];
        for (int j = 1; j < K; ++j)
        {
            if (mat_distance_fragment[i * nb_nodes + medoids[j]] < min)
            {
                min = mat_distance_fragment[i * nb_nodes + medoids[j]];
            }
        }
        #pragma omp atomic 
        sum += min;
    }
    return sum;
}
```

## MAIS

lors de test je me suis apercu que cette parallélisation entraine de trop grandes pertes de performances (2.5 fois pire). 

- Avec : temps total : 88286.6ms
- Sans : temps total : 34721.5ms

# `scatteredFloydAlgorithm` : boucle principales

La double boucle principale (celle qui fait les calculs de minimum pour chaque case du bloc) peut être parallélisé car chaque case est indépendantes des autres.

```c++
[...]

#pragma omp parallel for shared(bloc, recv_from_COMM_COL, recv_from_COMM_LINE)
        for (int i = 0; i < b; ++i)
        {
            for (int j = 0; j < b; ++j)
            {
                bloc[i * b + j] = min(bloc[i * b + j], recv_from_COMM_COL[j] + recv_from_COMM_LINE[i]);
            }
        }

[...]
```