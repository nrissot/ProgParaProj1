#ifndef FUNCP_HPP
#define FUNCP_HPP

#include <vector>
#include <math.h>
#include <ctime>
#include <cstdlib>

#include <mpi.h>
#include "ForGraph.hpp"

using namespace std;

/**
 * @brief permet de mettre en forme la matrice d'adjacence pour simplifier le MPI_Scatter. Cette fonction 
 * met aussi les 0 qui ne sont pas sur les diagonales à INF.
 * 
 * @param nb_nodes le nombre de noeuds du graphe.
 * @param mat_adjacence la matrice d'adjacence (de taille nb_nodes × nb_nodes).
 * @return int* une copie de la matrice d'adjacence préparée pour le scatter.
 */
int* prepareForScatter(int nb_nodes, int* mat_adjacence);

/**
 * @brief permet de remettre en forme la matrice reçue pour la reconvertir en matrice d'adjacence 
 * (opération inverse de prepareForScatter).
 * 
 * @param nb_nodes le nombre de noeuds du graphe.
 * @param mat_adjacence la matrice reçue (de taille nb_nodes × nb_nodes).
 * @return int* la matrice d'adjacence remise en forme.
 */
int* repareAfterGather(int nb_nodes, int* gathered_mat);

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


/**
 * @brief Permet de calculer les k meilleurs candidats de médoïdes locaux (ie sur un fragment de la matrice) 
 * pour accelerer l'éxecution de PAM. en ne choississant pas des candidats aléatoires.
 * 
 * @param mat_distance_fragment un fragment de la matrice des distances découpée en ligne (de taille nb_lignes_fragment × nb_nodes).
 * @param k le nombre de médoïde candidats a retourner (k <= nb_nodes).
 * @param nb_nodes le nombre total de noeud dans la matrice (taille des lignes de la matrice).
 * @param nb_lignes_fragment le nombre de noeud dont on a reçu les distances (nombre de lignes du fragment).
 *
 * @return un tableau d'entier avec pour valeur 0 ou 1. 
 * L'index i vaut 1 (resp. 0) si le i-ème noeud a été choisis (resp. n'a pas été choisit) comme médoïde.
 * Le tableau retourné aura exactement k valeurs 1 et nb_nodes - k valeurs 0
 */
int* findLocalMedoidCandidate(int* mat_distance_fragment, int k, int nb_nodes, int nb_lignes_fragment);

/**
* @brief  Traiter le tableau reduced des candidats locaux.
* @param  dadastruct un vecteur temporaire pour le traitement.
* @param  data_to_process le tableau réduit des candidats locaux.
* @param  nb_medoides le nombre de medoides voulu.
* @param  nb_noeud le nombre de noeud de notre matrice de distance.
* @return tableau des candidats globaux choisis.
 */
int* process_candidates(std::vector<std::vector<int>>* dadastruct,int* data_to_process,int nb_noeud,int nb_medoides);

/**
* @brief  Calcul le cout du swap d'un ménoide par un autre noeud du graphe.
* @param  current_candidates La liste des menoide a l'état courrant.
* @param  cout_locaux le tableeau contenant les couts de chaque swap possible.
* @param  data_fragment un fragment de la matrice des distances découpée en ligne (de taille nb_lignes_fragment × nb_nodes).
* @param  k le nombre de medoides voulu.
*/
void calcul_cout_swap(int* current_candidates, int* cout_locaux,int* data_fragment, int k);  // TODO: Just do it.

/**
* @brief  Calcule les meilleurs candidat sur une itération
* @param  cout_global_reduced La somme de tous les couts de swap locaux.
* @param  k le nombre de medoides voulu.
* @param  current_candidates La liste des menoide a l'état courrant.
* @return 0 si les menoid ne changent pas. 1 si il y a un changement.
*/
int choix_nouveaux_candidats(int* cout_global_reduced, int k, int* current_candidates); // TODO: Just do it.
#endif //FUNCP_HPP