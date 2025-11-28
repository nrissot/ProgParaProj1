#include <iostream>
#include <string.h>

#include <mpi.h>

#include "ForARN.hpp"
#include "Utils.hpp"
#include "FuncP.hpp"

using namespace std;

int main(int argc, char* argv[]) {
	if (argc != 4) {
		cout << "Usage : ./arn <fichier.fa> <k> <nb_seq>" << endl;
		return EXIT_FAILURE;
	}

	const int K = atoi(argv[2]);
	const int NB_SEQ = atoi(argv[3]);

	int pid, nprocs;
	int root = 0;

	// au cas ou on voudrait gerer des séquences dupliquées.
	int nb_nodes = NB_SEQ;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	char* sequences;

	if (pid == root) {
        sequences = readArnFromFile(argv[1], NB_SEQ);
	}

    
	int *displs, *sendcount;
	int nb_lignes = (pid < nb_nodes%nprocs) ? ((nb_nodes/nprocs)+1) : (nb_nodes/nprocs);
	int recvcount = nb_lignes * TAILLESEQ; 
	char* sequences_fragment = new char[recvcount];
    
	if (pid == root) {
		displs = new int[nprocs];
		sendcount = new int[nprocs];
		int displacement_index = 0;
		int nprocs_overload = nb_nodes % nprocs;

		for (int i = 0; i < nprocs; ++i) {
			int current_send_count_value = (i < nprocs_overload) ? (((nb_nodes/nprocs)+1)*TAILLESEQ) : ((nb_nodes/nprocs)*TAILLESEQ);
			sendcount[i] = current_send_count_value;
			displs[i] = displacement_index;
			displacement_index += current_send_count_value;
		}
	}
    
	MPI_Scatterv(sequences, sendcount, displs, MPI_CHAR, sequences_fragment, recvcount, MPI_CHAR, root, MPI_COMM_WORLD);
    
	if (pid == root) {
        delete[] displs;
		delete[] sendcount;
	}

    
    
	// déclaration de la matrice d'adjacence. 
	// Chaque Proc va avoir un tableau complet qui sera agrégé par root avec un reduce sum.
	// on initialise la matrice à 0 pour ne pas avoir à set les diagonales et les absences de liens.
	int* matrice_adjacence = new int[nb_nodes*nb_nodes]{0};

	// we add an extra line that may not be used for the communications with procs that 
	// dont have an overload.
	char* recv_buffer = new char[((nb_nodes/nprocs)+1) * TAILLESEQ];

    MPI_Barrier(MPI_COMM_WORLD);
    
	int next_nb_lignes;
	int next_recvcount;
    int next_offset;
    int self_offset = (pid < nb_nodes%nprocs) ? (pid * ((nb_nodes/nprocs)+1)) : ((nb_nodes%nprocs) * ((nb_nodes/nprocs)+1) + (pid - (nb_nodes%nprocs)) * (nb_nodes/nprocs));

    // ( i  <= nb_nodes%nprocs) ? ( i  * ((nb_nodes/nprocs)+1)) : ((nb_nodes%nprocs) * ((nb_nodes/nprocs)+1) + ( i  - (nb_nodes%nprocs)) * (nb_nodes/nprocs));
    // (pid <= nb_nodes%nprocs) ? (pid * ((nb_nodes/nprocs)+1)) : ((nb_nodes%nprocs) * ((nb_nodes/nprocs)+1) + (pid - (nb_nodes%nprocs)) * (nb_nodes/nprocs));
	for (int i = 0; i < nprocs; ++i) {
        if (pid == i) {
            // copy the fragment to the broadcasted buffer
			memcpy(recv_buffer, sequences_fragment, recvcount*sizeof(char));
			// fill in the amount broadcasted
			next_nb_lignes = nb_lignes;
			next_recvcount = recvcount;
            next_offset = self_offset;
		} else {
            // calculate the amount of data that will be transferred
			// ie, nb_lignes and recvcount value for proc n°i
			next_nb_lignes = (i < nb_nodes%nprocs) ? ((nb_nodes/nprocs)+1) : (nb_nodes/nprocs);
			next_recvcount = next_nb_lignes * TAILLESEQ;
            next_offset = (i < nb_nodes%nprocs) ? (i * ((nb_nodes/nprocs)+1)) : ((nb_nodes%nprocs) * ((nb_nodes/nprocs)+1) + (i - (nb_nodes%nprocs)) * (nb_nodes/nprocs));
		}
        
        MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(recv_buffer, next_recvcount, MPI_CHAR, i, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

		// if (pid != i) {
            // ajout des distance à la matrice d'adjacence
        buildMatrixFragment(matrice_adjacence, nb_nodes, sequences_fragment,self_offset, nb_lignes, recv_buffer, next_offset, next_nb_lignes);
		// }
		// comme la matrice est initialisée à 0, et que buildMatrixFragment prend en charge 
		// les ∞, il n'est pas nécessaire de remettre des valeurs pour la diagonale.
		// -> pas de else
		// --> pid n°i ne travaille pas, temps perdu ?
	} 


	// Agrégation des fragments de matrices construits.
	int *temp;
	if (pid == root) {
		temp = new int[nb_nodes*nb_nodes];
	}

	MPI_Reduce(matrice_adjacence, temp, nb_nodes*nb_nodes, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
    
    // if (pid == root) {
    //     // TODO:DELETEME:DEBUG
    //     cout << "matrice d'adjacence" << endl;
    //     affichage(temp, nb_nodes, nb_nodes, 3, INF);
    // }

    int* matrice_prep;

	if (pid == root) {
		matrice_prep = prepareForScatter(nb_nodes, temp);
		delete[] temp;
	}
	delete[] recv_buffer;



	// Algorithme de Roy-Floyd-Warshall pour calculer la matrice des distances
    int nb_blocs_par_lignes = (int) sqrt(nprocs);
    int b = nb_nodes / nb_blocs_par_lignes;

    if (pid == root) {
        if ((nb_nodes % b) != 0) {
            cout << "error: invalid value, (nb_nodes % b) != 0" << endl;
            exit(1);
        }
        if ((nb_nodes/b * nb_nodes/b) != (nb_blocs_par_lignes * nb_blocs_par_lignes)) {
            cout << "error: invalid value, (nb_nodes/b * nb_nodes/b) != (nb_blocs_par_lignes * nb_blocs_par_lignes)" << endl; 
            exit(1);
        }
    }

    int *bloc = new int[b*b]();
    
    MPI_Scatter(matrice_prep, b*b, MPI_INT, bloc, b*b, MPI_INT, root, MPI_COMM_WORLD);

    // Creation des Communicateur ligne et colonne
    int* ndims = new int[2]{nb_blocs_par_lignes, nb_blocs_par_lignes};
    int* periods = new int[2]{0,0};

    // Definitions des dimensions à garder pour la création des communicateurs par Cart_sub
    int* remain_dims_col = new int[2]{0,1};
    int* remain_dims_line = new int[2]{1,0};

    MPI_Comm MPI_COMM_CART;
    MPI_Comm MPI_COMM_COL;
    MPI_Comm MPI_COMM_LINE;

    
    MPI_Cart_create(MPI_COMM_WORLD, 2, ndims, periods, false, &MPI_COMM_CART);

    MPI_Cart_sub(MPI_COMM_CART, remain_dims_col, &MPI_COMM_COL);
    MPI_Cart_sub(MPI_COMM_CART, remain_dims_line, &MPI_COMM_LINE);         


    // Libération des divers array utilisés pour la séparation du 
    // communicateur world en grille.
    delete[] ndims;
    delete[] periods;
    delete[] remain_dims_line;
    delete[] remain_dims_col;

    scatteredFloydAlgorithm(bloc, b, nb_nodes, MPI_COMM_COL, MPI_COMM_LINE);

    int *mat_gathered, *mat_distances;

    if (pid == root) {
       mat_gathered = new int[nb_nodes*nb_nodes];
    }

    MPI_Gather(bloc, b*b, MPI_INT, mat_gathered, b*b, MPI_INT, root, MPI_COMM_WORLD);
    delete[] bloc;
    
    if (pid == root) { 
        mat_distances = repareAfterGather(nb_nodes, mat_gathered);
    }

    // if (pid == root) {
    //     // TODO:DELETEME:DEBUG
    //     cout << "\n\n\nmatrice distances" << endl;
    //     affichage(mat_distances, nb_nodes, nb_nodes, 3, INF);
    // }


	// PAM pour calculer les k-médoïdes
    
    nb_lignes = (nb_nodes%nprocs >= pid) ? ((nb_nodes/nprocs)+1) : (nb_nodes/nprocs);
    recvcount = nb_lignes * nb_nodes; 
    int* mat_distances_fragment = new int[recvcount];
    
    if (pid == root) {
        displs = new int[nprocs];
        sendcount = new int[nprocs];
        int displacement_index = 0;
        int nprocs_overload = nb_nodes % nprocs;

        for (int i = 0; i < nprocs; ++i) {
            int current_send_count_value = (nprocs_overload > i) ? (((nb_nodes/nprocs)+1)*nb_nodes) : ((nb_nodes/nprocs)*nb_nodes);
            sendcount[i] = current_send_count_value;
            displs[i] = displacement_index;
            displacement_index += current_send_count_value;
        }
    }

    
    MPI_Scatterv(mat_distances, sendcount, displs, MPI_INT, mat_distances_fragment, recvcount, MPI_INT, root, MPI_COMM_WORLD);
    
    if (pid == root) {
        delete[] displs;
        delete[] sendcount;
    }
    
    int* local_chosen_candidates = findLocalMedoidCandidate(mat_distances_fragment, K, nb_nodes, (nb_nodes / nprocs));
    
    int* reduced_candidates;
    
    if (pid == root) {
        reduced_candidates = new int[nb_nodes];
    }
    
    MPI_Reduce(local_chosen_candidates,reduced_candidates,nb_nodes,MPI_INT,MPI_SUM,root,MPI_COMM_WORLD);

    // if (pid == root) {
    //     // TODO:DELETEME:DEBUG
    //     cout << "\n\n\nlocal chosen candidates (reduced)" << endl;
    //     affichage(reduced_candidates, 1,nb_nodes, 3, INF);
    // }
    
	delete[] local_chosen_candidates;
    
    int *medoids = new int[K];
    
    if (pid == root) {
        //                    max amount of votes a medoid can receive ↓ (if it shoot to high, can add a max() pass before).
        medoids = get_k_best_elt(reduced_candidates, nb_nodes, K, nprocs-1);
    }
    
    MPI_Bcast(medoids, K, MPI_INT, root, MPI_COMM_WORLD);
    
    int start_cost = calculate_cost_fragment(medoids, K, mat_distances_fragment, nb_nodes, nb_lignes);
    int global_cost = 0;
    
    MPI_Reduce(&start_cost, &global_cost, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
    
    
    // if (pid == root) {
    //     // TODO:DELETEME:DEBUG
    //     cout << "starting cost : " << global_cost << endl ;
    //     affichage(medoids, 1, K, 3, INF);
    // }
    
    int *permutations = new int[K*K];
    int *permutation_cost = new int[K];
    temp = new int[K];


    for (int i=0; i < nb_nodes; ++i) {
        if (pid == root) {
            // generate the permutations
            if (!is_in(medoids, i, K)) {
                for (int j = 0; j < K; ++j) {
                    for (int k=0; k < K; ++k) {
                        permutations[j*K + k] = ((j == k) ? i : medoids[k]);
                    }
                }
            }
        }
        MPI_Bcast(permutations, K*K, MPI_INT, root, MPI_COMM_WORLD);

        for (int j = 0; j < K; ++j) {
            permutation_cost[j] = calculate_cost_fragment(&permutations[j*K], K, mat_distances_fragment, nb_nodes, nb_lignes);
        }


        MPI_Reduce(permutation_cost, temp, K, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);

        memcpy(permutation_cost, temp, K * sizeof(int));
        

        if (pid == root) {
            int best_cost_index = min_elt_index(permutation_cost, K);
            int new_cost = permutation_cost[best_cost_index];
            if (new_cost < global_cost) {
                // make the permutation
                global_cost = new_cost;
                memcpy(medoids, &permutations[best_cost_index*K], K * sizeof(int));
            }
        }

        // TODO:OPTIMIZE:MAYBE: can this be optimized ? how can we know when we dont need to bcast ?
        MPI_Bcast(medoids, K, MPI_INT, root, MPI_COMM_WORLD);
    }

    if (pid == root) {
        // TODO:DELETEME:DEBUG
        cout << "ending cost : " << global_cost << endl ;
        affichage(medoids, 1, K, 3, INF);
    }

	// Export (+ calcul des communautées)

    // final cleanup
    delete[] medoids;
    delete[] permutations;
    delete[] permutation_cost;
    delete[] temp;
    delete[] mat_distances_fragment;
    
    if (pid == root) {
        delete[] mat_gathered;
        delete[] mat_distances;
        // delete[] reduced_candidates;
		delete[] sequences;
	}
	delete[] sequences_fragment;
    delete[] matrice_adjacence;

	MPI_Finalize();
	return 0;
}