#include "FuncP.hpp"
#include <iostream>

int* prepareForScatter(int nb_nodes,int *mat_adjacence) {
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    
    int nb_blocs_par_lignes = int(sqrt(double(nprocs)));
    int b = nb_nodes / nb_blocs_par_lignes;

    // matrice renvoyée.
    int* resultat = new int[nb_nodes*nb_nodes]();

    // position de la "tête d'écriture" dans résultat.
    int indice = 0;

    // On parcours les blocs, (chaque processeur reçevra un blocs, nb_blocs = nprocs).
    for (int bloc_id = 0; bloc_id < nprocs; ++bloc_id) {
        int numero_colonne_bloc = bloc_id % nb_blocs_par_lignes;
        int numero_ligne_bloc = bloc_id / nb_blocs_par_lignes;

        // flag permettant de savoir si notre bloc contient la diagonal de la matrice d'origine.
        int bloc_diagonal = numero_colonne_bloc == numero_ligne_bloc ? 1 : 0;
        // position du premier element du bloc dans la matrice d'adjacence.
        int position_debut_bloc = (numero_ligne_bloc * (b*b) * nb_blocs_par_lignes) + (numero_colonne_bloc * b);

        // Pour chaque element du bloc
        for (int x = 0; x < b; ++x) {
            for (int y = 0; y < b; ++y) {
                // On lit l'element a ecrire dans le résultat
                // + (x * nb_nodes) permet de se décaler d'une ligne vers le bas dans la matrice d'adjacence.
                int value = mat_adjacence[position_debut_bloc + (x * nb_nodes) + y];

                // Logique de remplacement de la valeur 0 par l'infini lorsqu'il n'y a pas de lien entre les noeuds.
                if(value == 0){
                    // Si x == y dans un bloc contenant la diagonale globale on met à 0 
                    // (indique qu'un iteration est à distance 0 de lui même).
                    if(bloc_diagonal && x==y) {
                        value = 0;
                    } else {
                        value = INF;
                    } 
                }
            resultat[indice] = value;
            ++indice;
            }
        }
    }
    return resultat;
}

int *repareAfterGather(int nb_nodes, int *gathered_mat) {
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    
    int nb_blocs_par_lignes = int(sqrt(double(nprocs)));
    int b = nb_nodes / nb_blocs_par_lignes;

    int* resultat = new int[nb_nodes*nb_nodes];

    int reading_index = 0;
    // pour chaque bloc
    for (int block_id = 0; block_id < nprocs; ++block_id) {
        int x = (block_id / nb_blocs_par_lignes);
        int y = (block_id % nb_blocs_par_lignes);
        int writing_index = (x * nb_nodes * b) + (y * b);
        // pour chaque indice dans le bloc
        for (int i = 0; i < b; ++i) {
            for (int j = 0; j < b; ++j) {
                resultat[writing_index] = gathered_mat[reading_index];
                reading_index++;
                writing_index++;
            }
            writing_index += (nb_blocs_par_lignes -1)*b;
        }
    }
    return resultat;
}

void scatteredFloydAlgorithm(int* bloc, int b, int nb_nodes, MPI_Comm MPI_COMM_COL, MPI_Comm MPI_COMM_LINE) {
    int col_pid, line_pid, pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_rank(MPI_COMM_COL, &col_pid);
    MPI_Comm_rank(MPI_COMM_LINE, &line_pid);

    int* recv_from_COMM_COL = new int[b];   // taille (b*1)
    int* recv_from_COMM_LINE = new int[b];  // taille (1*b)

    for (int l=0; l < nb_nodes; ++l) {
        int position_pid_dans_COMM = l / b;
        int position_dans_bloc = l % b;

        if (col_pid == position_pid_dans_COMM) {
            // on copie les element de la ligne n° position du bloc dans recv from COL
            for (int i = 0; i < b; ++i) {
                recv_from_COMM_COL[i] = bloc[i*b + position_dans_bloc];
            }
        }
        if (line_pid == position_pid_dans_COMM) {
            // on copie les element de la colonne n° position du bloc dans recv from LINE
            for (int i = 0; i < b; ++i) {
                recv_from_COMM_LINE[i] = bloc[position_dans_bloc*b + i];
            }
        }
        MPI_Bcast(recv_from_COMM_COL, b, MPI_INT, position_pid_dans_COMM, MPI_COMM_COL);
        MPI_Bcast(recv_from_COMM_LINE, b, MPI_INT, position_pid_dans_COMM, MPI_COMM_LINE);

        for (int i = 0; i < b; ++i) {
            for (int j = 0; j < b; ++j) {
                bloc[i*b +j] = min(bloc[i*b +j], recv_from_COMM_LINE[j]+recv_from_COMM_COL[i]);
            }
        }

    }
    delete[] recv_from_COMM_COL;
    delete[] recv_from_COMM_LINE;
}
