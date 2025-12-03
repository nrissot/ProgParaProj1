#include <iostream>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <chrono>

#include "ForARN.hpp"
#include "Utils.hpp"
#include "FuncP.hpp"

using namespace std;

int main(int argc, char* argv[]) {
	if (argc != 4) {
		cout << "Usage : ./arn <fichier.fa> <k> <nb_seq>" << endl;
		return EXIT_FAILURE;
	}

    auto t_start = std::chrono::high_resolution_clock::now();

	const int K = atoi(argv[2]);
	const int NB_SEQ = atoi(argv[3]);

	int pid, nprocs;
	int root = 0;

    // used for the exported dotfile community detection creation
    string colors[4] = {"aquamarine", "goldenrod1", "hotpink1", "chartreuse"};

	// au cas ou on voudrait gerer des séquences dupliquées.
	int nb_nodes = NB_SEQ;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    auto t_init = std::chrono::high_resolution_clock::now();

	char* sequences;

	if (pid == root) {
        sequences = readArnFromFile(argv[1], NB_SEQ);
	}
    auto t_lecture = std::chrono::high_resolution_clock::now();

    
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

    
	int next_nb_lignes;
	int next_recvcount;
    int next_offset;
    int self_offset = (pid < nb_nodes%nprocs) ? (pid * ((nb_nodes/nprocs)+1)) : ((nb_nodes%nprocs) * ((nb_nodes/nprocs)+1) + (pid - (nb_nodes%nprocs)) * (nb_nodes/nprocs));

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
        
		MPI_Bcast(recv_buffer, next_recvcount, MPI_CHAR, i, MPI_COMM_WORLD);

        buildMatrixFragment(matrice_adjacence, nb_nodes, sequences_fragment,self_offset, nb_lignes, recv_buffer, next_offset, next_nb_lignes);
	} 


	// Agrégation des fragments de matrices construits.
	int *temp;
	if (pid == root) {
		temp = new int[nb_nodes*nb_nodes];
	}

	MPI_Reduce(matrice_adjacence, temp, nb_nodes*nb_nodes, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
    auto t_matrice_adjacence = std::chrono::high_resolution_clock::now();

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
    int* remain_dims_col = new int[2]{1,0};
    int* remain_dims_line = new int[2]{0,1};

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
    delete[] remain_dims_col;
    delete[] remain_dims_line;

    scatteredFloydAlgorithm(bloc, b, nb_nodes, MPI_COMM_LINE, MPI_COMM_COL);

    int *mat_gathered, *mat_distances;

    if (pid == root) {
       mat_gathered = new int[nb_nodes*nb_nodes];
    }

    MPI_Gather(bloc, b*b, MPI_INT, mat_gathered, b*b, MPI_INT, root, MPI_COMM_WORLD);
    delete[] bloc;
    
    if (pid == root) { 
        mat_distances = repareAfterGather(nb_nodes, mat_gathered);
    }

    auto t_fin_floyd = std::chrono::high_resolution_clock::now();

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
    
    int* local_chosen_candidates;

    if (nprocs > 1) {
        local_chosen_candidates = findLocalMedoidCandidate(mat_distances_fragment, K, nb_nodes, (nb_nodes / nprocs));
    } else  {
        // special test case for a sequential execution (mpirun -np 1) of the algorithm
        local_chosen_candidates = new int[K];
        srand(time({})); // use current time as seed for the RNG.
        for (int i = 0; i < K; ++i) {
            int choice = rand() % nb_nodes;
            int flag = 0;
            // check if the chosen medoid was not already chosen.
            for (int j = 0; j < i; ++j) {
                if (local_chosen_candidates[j] == choice) {
                    flag =1;
                    break;
                }
            }
            if (flag == 0) {
                local_chosen_candidates[i] = choice;
            } else {
                // we have to try again
                --i;
            }
        }
    }
    int* reduced_candidates;
    
    if (pid == root) {
        reduced_candidates = new int[nb_nodes];
    }
    
    MPI_Reduce(local_chosen_candidates,reduced_candidates,nb_nodes,MPI_INT,MPI_SUM,root,MPI_COMM_WORLD);
    
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

        if (pid == root) {
            memcpy(permutation_cost, temp, K * sizeof(int));
            int best_cost_index = min_elt_index(permutation_cost, K);
            int new_cost = permutation_cost[best_cost_index];
            if (new_cost < global_cost) {
                // make the permutation
                global_cost = new_cost;
                memcpy(medoids, &permutations[best_cost_index*K], K * sizeof(int));
            }
        }

        MPI_Bcast(medoids, K, MPI_INT, root, MPI_COMM_WORLD);
    }
    auto t_fin_pam = std::chrono::high_resolution_clock::now();

    if (pid == root) {
        cout << "Coût final : " << global_cost << "\nMédoïdes : " << endl ;
        affichage(medoids, 1, K, 3, INF);
        cout << "---" << endl;
    }

	// Export (+ calcul des communautées) 
    // NODES
    // C2N20 [label="C2N20",color="darkorchid"];
    if (pid == root) {
        ofstream ofile;
        ofile.open("output/arn.dot");
        ofile.clear();
        ofile << "graph output_arn {\nnode [shape=circle, style=filled, color=lightyellow, fontcolor=black];\nedge [color=black, fontcolor=blue];\n" << endl;
        
        for (int i = 0; i < nb_nodes; ++i) {
            // find the closest medoid (ie, community)
            int closest_medoid_idx = 0;
            int closest_medoid_dist = mat_distances[i*nb_nodes + medoids[0]];
            for (int k = 1; k < K; ++k) {
                if (mat_distances[i*nb_nodes + medoids[k]] < closest_medoid_dist) {
                    closest_medoid_idx = k;
                    closest_medoid_dist = mat_distances[i*nb_nodes + medoids[k]];
                }
            }
            //                                                    reused colors if k > 4 ↓
            ofile << i << " [label=\"" << i << "\",color=\"" << colors[closest_medoid_idx%4] ;

            // if the node is a medoid
            if (closest_medoid_dist == 0) {
                // make the node bigger and a pentagon instead of an ellipse
                ofile << "\",fixedsize=true,width=2,height=2,shape=\"pentagon";
            }

            ofile << "\"];" << endl; ;

            // find the edges
            // C2N21 -- C2N13 [label="1", weight=1];
            for (int j = i; j < nb_nodes; ++j) {
                const int w = mat_distances[i*nb_nodes + j];
                if (w <= 70 && w != 0) {
                    ofile << i << " -- " << j << " [label=\"" << w << "\",weight=" << w << "];" << endl ;
                }
            }
        }
        ofile << "}" << endl;
        ofile.close();
    }

    // final cleanup
    delete[] medoids;
    delete[] permutations;
    delete[] permutation_cost;
    delete[] temp;
    delete[] mat_distances_fragment;
    
    if (pid == root) {
        delete[] mat_gathered;
        delete[] mat_distances;
		delete[] sequences;
	}
	delete[] sequences_fragment;
    delete[] matrice_adjacence;

	MPI_Finalize();
    if(pid == root)
    {
        using std::chrono::duration;
        using std::chrono::duration_cast;
        using std::chrono::high_resolution_clock;
        using std::chrono::milliseconds;

        duration<double, std::milli> ms_double_total = t_fin_pam - t_start;
        std::cout<<"temps total : " << ms_double_total.count() << "ms" << endl;
        
        duration<double, std::milli> ms_double_init = t_init - t_start;
        std::cout<<"temps init : " << ms_double_init.count() << "ms" << endl;

        duration<double, std::milli> ms_double_lecture = t_lecture - t_init;
        std::cout<<"temps lecture : " << ms_double_lecture.count() << "ms" << endl;

        duration<double, std::milli> ms_double_distance = t_matrice_adjacence - t_lecture;
        std::cout<<"temps matrice distance : " << ms_double_distance.count() << "ms" << endl;

        duration<double, std::milli> ms_double_floyd = t_fin_floyd - t_matrice_adjacence;
        std::cout<<"temps floyd : " << ms_double_floyd.count() << "ms" << endl;

        duration<double, std::milli> ms_double_pam = t_fin_pam - t_fin_floyd;
        std::cout<<"temps pam : " << ms_double_pam.count() << "ms" << endl;
    }
	return 0;
}