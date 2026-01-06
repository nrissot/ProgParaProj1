#ifndef NEEDLEMAN
#define NEEDELMAN

#include <algorithm> // needed for max()
#include <omp.h>

using namespace std;

#define MATCH_REWARD 1
#define MISMATCH_PENALTY -1
#define GAP_PENALTY -3
#define EXTEND_PENALTY -1

/**
 * @brief Permet de calculer la distance de Needleman-Wunsch telle que décrite dans
 *        https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm de maniere séquentielle.
 *
 * @param A la premiere chaîne de caractères
 * @param lenA la taille de A
 * @param B la seconde chaîne de caractères
 * @param lenB la taille de B
 *
 * @return le score (en int) de la comparaison entre A et B
 */
int nw_distance_seq(const char *A, int lenA, const char *B, int lenB);

/**
 * @brief Permet de calculer la distance de Needleman-Wunsch telle que décrite dans
 *        https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm de maniere
 *        parallélisée en utilisant OpenMP
 *
 * @param A la premiere chaîne de caractères
 * @param lenA la taille de A
 * @param B la seconde chaîne de caractères
 * @param lenB la taille de B
 *
 * @return le score (en int) de la comparaison entre A et B
 */
int nw_distance_omp(const char *A, int lenA, const char *B, int lenB);

/**
 * @brief Fonction Helper permettant de réutiliser le comportement d'écriture et de test dans les méthodes
 *        séquentielles et parallélisée. Teste laquelle des trois options (CF Needleman/Needleman–Wunsch algorithm - Wikipedia.pdf)
 *        est la plus adaptée, et ecrit le score calculé en conséquence dans F, et si le score résulte d'un trou dans H
 *
 * @param F la matrice des scores
 * @param H la matrice des trou
 * @param i l'indice i courant
 * @param j l'indice j courant
 * @param a le caractère de A qui est comparé (A[i])
 * @param b le caractère de B qui est comparé (B[j])
 * @param lenA la longueur de la chaine A, ie la hauteur des matrices F et H
 * @param lenB la longueur de la chaine B, ie la largeur des matrices F et H
 */
void _choose_and_write(int *F, bool *H, int i, int j, char a, char b, int lenA, int lenB);

#endif