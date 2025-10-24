#ifndef FUNCP_HPP
#define FUNCP_HPP

#include <vector>
#include <math.h>

using namespace std;

/**
 * @brief 
 * 
 * @param nb_nodes le nombre de noeuds du graphe.
 * @param mat_adjacence la matrice d'adjacence (de taille nb_nodes × nb_nodes).
 * @param nprocs le nombre de processeurs.
 * @return int* une copie de la matrice d'adjacence préparée pour le scatter.
 */
int* prepareForScatter(int nb_nodes, int* mat_adjacence, int nprocs, int inf);

#endif //FUNCP_HPP