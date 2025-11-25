#include "FuncP.hpp"
#include <iostream>  // TODO:DELETEME
#include "Utils.hpp" // TODO:DELETEME

// ____ PRIVATE HELPERS ________________________________________________________________________________________________

/**
 * @brief fonction helper privée pour verifier la présence d'un element dans un tableau.
 * 
 * @param tab le tableau (de taille len).
 * @param elt l'élément recherché.
 * @param len la taille de tab.
 * 
 * @return 1 si l'élément est présent, 0 sinon.
 */
int is_in(int* tab, int elt, int len) {
    for (int i = 0; i < len; ++i) {
        if (tab[i] == elt) {
            return 1;
        }
    }
    return 0;
}

/**
 * @brief fonction helper privée qui permet de calculer la somme des poids des k candidats choisit.
 * 
 * @param mat_distance_fragment un fragment de la matrice des distances découpée en ligne (de taille nb_lignes_fragment × nb_nodes).
 * @param candidates les k candidats choisi
 * @param k le nombre de médoïde candidats a retourner (k <= nb_nodes).
 * @param nb_nodes le nombre total de noeud dans la matrice (taille des lignes de la matrice).
 * @param nb_lignes_fragment le nombre de noeud dont on a reçu les distances (nombre de lignes du fragment).
 *
 * 
 * @return la somme des poids des k candidats choisit.
 */
int cost_from_candidate_set(int *mat_distance_fragment, int* candidates, int k, int nb_nodes, int nb_lignes_fragment) {
    int sum = 0;
    for(int current_node = 0; current_node < nb_lignes_fragment; current_node++) {
        int min_current = INF;          
        for (int i = 0 ; i < k; i++) {
            // NOTE: old way to spot error, shoudlnt happen TODO:DELETEME ?
            // if (mat_distance_fragment[current_node*nb_nodes + candidates[i]] == 0) {
            //     return -1;
            // }
            min_current = min(mat_distance_fragment[current_node*nb_nodes + candidates[i]], min_current);
        }
        sum += min_current;
    }
    return sum;
}

// ____ HPP FUNCTIONS _________________________________________________________________________________________________

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

int *findLocalMedoidCandidate(int *mat_distance_fragment, int k, int nb_nodes, int nb_lignes_fragment) {
    int* candidates = new int[k]{};
    srand(time({})); // use current time as seed for the RNG.
    
    // select k unique candidates.
    for (int i = 0; i < k;) {
        int used = 0;
        int choice = rand() % nb_nodes;
        for(int current_line = 0; current_line<nb_lignes_fragment;current_line++ ){
            if(mat_distance_fragment[current_line*nb_nodes+choice] == 0){
                used = 1;break;
            }
        }
        for (int j = 0; j < i; ++j) {
            if (candidates[j] == choice){
                used = 1;
                break;
            }
        }
        if (used == 0) {
            candidates[i] = choice;
            ++i;
        }
    }

    // find best local medoid.
    int flag;
    do {
        flag = 0;
        // current cost of the k chosen medoid.
        int current_cost = cost_from_candidate_set(mat_distance_fragment,candidates,k,nb_nodes,nb_lignes_fragment);
        
        int* copy = new int[k];
        // Array of local index node in our matrix fragment.
        int* index_locaux = new int[nb_lignes_fragment];

        
        for (int i = 0; i < nb_lignes_fragment; ++i) {
            for (int j = 0; j < nb_nodes; ++j) {
                if(mat_distance_fragment[i*nb_nodes + j] == 0) {
                    // if the value is the distance from a node to itself, (ie, value=0)
                    // we add it as a local node.
                    index_locaux[i] = j;
                    break;
                }
            }
        }
        int pid; 
        MPI_Comm_rank(MPI_COMM_WORLD, &pid);

        // if(pid == 15){
        //     cout << "affichage du fragment de matrice de 15: " << endl;
        //     affichage(mat_distance_fragment, nb_lignes_fragment, nb_nodes,2,1000);
        //     cout << "hmmmm trés intéressant hmm hmm " << endl; // TODO:DELETEME

        // }
        
        
        // for each medoid
        for (int i = 0; i < k; ++i) {
            
            // we try every node that is not already a medoid to see if we can lower the cost
            for (int j = 0; j < nb_nodes; ++j) {
                memcpy(copy,candidates,sizeof(int)*k);
                // rajouter "au bout" index_locaux

                if(!is_in(copy,j,k) && !is_in(index_locaux,j,nb_lignes_fragment)) { // && !is_in index_locaux
                    copy[i] = j;
                    int new_cost = cost_from_candidate_set(mat_distance_fragment,copy,k,nb_nodes,nb_lignes_fragment);
                    // if the candidate is not a member of the fragment
                    if(current_cost > new_cost){
                        // we can lower the cost, we update the candidates array.
                        candidates[i]=copy[i];

                        // we change something, we update the flag.
                        flag=1;                        
                    }
                }
            }
        }
        delete[] copy;
        
    } while (flag != 0);
    
    int* candidates_flags = new int[nb_nodes]{};
    //                                       ↑ initialise the values of the array to 0

    for (int i = 0; i < k; ++i) {
        candidates_flags[candidates[i]] = 1;
    }
    delete[] candidates;
    return candidates_flags;
}

int* process_candidates(std::vector<std::vector<int>>* datastruct,int* data_to_process,int nb_noeud,int nb_medoides){
    int medoides_added = 0; 
    int* data_processed = new int[nb_medoides];
    // on parcours data to process et on remplit datastruct.
    for(int i = 0;i<nb_noeud;++i){
        (*datastruct)[data_to_process[i]].push_back(i);
    }
    for (int i = (*datastruct).size() - 1; i >= 0; --i) {
        for(int j = 0 ; j < (*datastruct)[i].size();j++){
            data_processed[medoides_added] = (*datastruct)[i][j];
            medoides_added++;
            if(medoides_added==nb_medoides){
                break;
            }
        }
        if(medoides_added==nb_medoides){
            break;
        }
    }

    return data_processed;
}

void calcul_cout_swap(int *current_candidates, int *cout_locaux, int *mat_distance_fragment, int k, int nb_nodes, int nb_lignes_fragment) {
    int *new_candidates = new int[k];

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < nb_nodes; ++j) {
            // if the current index is not already a candidates
            if (!is_in(current_candidates, j, k)) {
                memcpy(new_candidates, current_candidates,k*sizeof(int));
                new_candidates[i] = j;
                cout_locaux[(i*nb_nodes)+j] = cost_from_candidate_set(mat_distance_fragment, new_candidates, k, nb_nodes, nb_lignes_fragment);
            }
        }
    }
}

int choix_nouveaux_candidats(int nb_nodes, int *cout_global_reduced, int k, int *current_candidates) {
    int flag = 0;
    for (int i = 0; i < k; ++i) {
        int optimal_cost = INF;
        for (int j = 0; j < nb_nodes; ++j) {
            if (!is_in(current_candidates,j,k) && cout_global_reduced[(i*nb_nodes)+j] < optimal_cost) {
                optimal_cost = cout_global_reduced[(i*nb_nodes)+j];
                current_candidates[i] = j;
                flag = 1;
            }
        }
    }
    return flag;
}
