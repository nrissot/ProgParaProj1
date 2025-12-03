#ifndef FORARN_HPP
#define FORARN_HPP

#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

/**
 * @brief taille des séquences d'ARN.
 */
#define TAILLESEQ 100

/**
 * @brief distance seuil au dela de laquelle on ne crée pas de lien entre deux séquences d'ARN.
 */
#define SEUILDISTANCE 70

/**
 * @brief ouvre puis parse un fichier .fa pour en générer le tableau contenant les chaines d'ARN.
 *
 * @param f le path du fichier .fa qui contient les séquences d'ARN
 * @param nb_nodes le nombre de séquences d'ARN contenues dans le tableau
 *
 * @return un tableau de caractères de taille nb_nodes * TAILLESEQ contenant les chaines d'ARN
 */
char* readArnFromFile(char* f, int nb_nodes);

/**
 * @brief permet de calculer la distance de hamming 2 à 2 entre les séquences d'ARN passées en argument et rentrée au bon endroit dans le fragment de la matrice d'adjacence 
 *
 * @param mat_adjacence [MODIFIÉ PAR LA FONCTION] matrice d'adjacence (de taille nb_nodes * nb_nodes).
 * @param nb_nodes le nombre total de séquences d'ARN.
 * @param sequence_fragment le fragment du tableau des séquence d'ARN reçu initialement par le processus (de taille nb_lignes * TAILLESEQ).
 * @param fragment_offset l'indice du premier élément de sequence_fragment dans la matrice d'adjacence.
 * @param nb_lignes le nombre de lignes contenues dans sequence_fragment.
 * @param recv_buffer le fragment du tableau des séquences d'ARN reçu par broadcast d'un autre proc (de taille recv_nb_lignes * TAILLESEQ).
 * @param recv_offset l'indice du premier élément de recv_buffer dans la matrice d'adjacence.
 * @param recv_nb_lignes le nombre de lignes contenues dans recv_buffer.
 */
void buildMatrixFragment(int* mat_adjacence, int nb_nodes, char* sequence_fragment, int fragment_offset, int nb_lignes, char* recv_buffer, int recv_offset, int recv_nb_lignes);
#endif