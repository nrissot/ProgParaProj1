#ifndef FUNCP_HPP
#define FUNCP_HPP

#include <vector>
#include <math.h>

#include <mpi.h>

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

/**
 * @brief effectue l'algorithme de Roy-Floyd-Warshall sur le fragment de la matrice d'adjacence dédié
 * 
 * @param bloc le fragment de la matrice d'adjacence correspondant au bloc (de taille b × b)
 * @param b la taille du bloc
 * @param nb_nodes le nombre total de noeuds du graphe
 * @param MPI_COMM_COL le communicateur en colonne
 * @param MPI_COMM_LINE le communicateur en ligne
 */
void scatteredFloydAlgorithm(int* bloc, int b, int nb_nodes, MPI_Comm MPI_COMM_COL, MPI_Comm MPI_COMM_LINE);

#endif //FUNCP_HPP