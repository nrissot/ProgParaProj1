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
int* prepareForScatter(int nb_nodes, int* mat_adjacence); //DEBUG:x

/**
 * @brief permet de remettre en forme la matrice reçue pour la reconvertir en matrice d'adjacence 
 * (opération inverse de prepareForScatter).
 * 
 * @param nb_nodes le nombre de noeuds du graphe.
 * @param mat_adjacence la matrice reçue (de taille nb_nodes × nb_nodes).
 * @return int* la matrice d'adjacence remise en forme.
 */
int* repareAfterGather(int nb_nodes, int* gathered_mat); //DEBUG:x

/**
 * @brief effectue l'algorithme de Roy-Floyd-Warshall sur le fragment de la matrice d'adjacence dédié
 * 
 * @param bloc le fragment de la matrice d'adjacence correspondant au bloc (de taille b × b)
 * @param b la taille du bloc
 * @param nb_nodes le nombre total de noeuds du graphe
 * @param MPI_COMM_COL le communicateur en colonne
 * @param MPI_COMM_LINE le communicateur en ligne
 */
void scatteredFloydAlgorithm(int* bloc, int b, int nb_nodes, MPI_Comm MPI_COMM_COL, MPI_Comm MPI_COMM_LINE); //DEBUG:x


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
int* findLocalMedoidCandidate(int* mat_distance_fragment, int k, int nb_nodes, int nb_lignes_fragment); //DEBUG:x

// NOUVELLES FONCTIONS NATHAN 2025.11.25
/**
 * @brief calcule et retourne les indices des k meilleurs élement du tableau passé en paramètre
 * @param tab le tableau (de taille length) des options
 * @param k le nombre d'élements voulus.
 * @param max la valeur maximum que peut prendre un élément
 * @return un tableau de taille k des meilleurs éléments.
 */
int* get_k_best_elt(int* tab, int length, int k, int max) ; //DEBUG:x

/**
 * @brief verifie la présence d'un element dans un tableau en temps O(len).
 * 
 * @param tab le tableau (de taille len).
 * @param elt l'élément recherché.
 * @param len la taille de tab.
 * 
 * @return true si l'élément est présent, false sinon.
 */
bool is_in(int* tab, int elt, int len); //DEBUG:x

/**
 * @brief retourne l'index du plus petit element du tableau en temps O(len).
 * 
 * @param tab le tableau (de taille len)
 * @param len la taille du tableau tab
 * 
 * @return l'index du plus petit element du tableau 
 */
int min_elt_index(int* tab, int len); //DEBUG:x

/**
 * @brief calcule le fragment de cout pour les médoïdes passé en paramètres en temps O(nb_lignes*k)
 * 
 * @param medoids le tableau des indices des médoïdes (de taille k)
 * @param k le nombre de médoides
 * @param mat_distance_fragment le fragment de la matrice des distances (de taille nb_nodes * nb_lignes)
 * @param nb_nodes le nombre de noeuds du graphe, la largeur de mat_distance_fragment
 * @param nb_lignes le nombre de lignes du fragement, la hauteur de mat_distance_fragment
 * 
 * @return la somme des distances des (nb_lignes) noeud du fragment avec leur plus proche médoïde 
 */
int calculate_cost_fragment(int* medoids, int k, int* mat_distance_fragment, int nb_nodes, int nb_lignes); //DEBUG:x

#endif //FUNCP_HPP