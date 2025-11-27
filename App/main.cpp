#include <iostream>
#include <string>
#include <map>
#include <vector>

#include <mpi.h>

#include "Utils.hpp"
#include "ForGraph.hpp"
#include "FuncP.hpp" 

using namespace std;


int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "Usage : ./main fichier.dot (graphe au format dot) k";
        return EXIT_FAILURE;
    }

    
    int pid, nprocs;
    int root = 0;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    
    int *mat_adjacence, *mat_preparee;
    
    int nb_nodes;
    
    
    if (pid == root) {
        char* file_name = argv[1];
        map<string,int> my_nodes; // Dans le .dot les sommets sont désignés par un nom 
        
        mat_adjacence = lectureGraphe(file_name,&nb_nodes,&my_nodes);        
        mat_preparee = prepareForScatter(nb_nodes, mat_adjacence);
    }
    
    MPI_Bcast(&nb_nodes, 1, MPI_INT, root, MPI_COMM_WORLD);

    int nb_blocs_par_lignes = (int) sqrt(nprocs);
    int b = nb_nodes / nb_blocs_par_lignes;

    if (pid == root) {
        if ((nb_nodes % b) != 0) {
            cout << "invalid value, (nb_nodes % b) != 0" << endl;
            exit(1);
        }
        if ((nb_nodes/b * nb_nodes/b) != (nb_blocs_par_lignes * nb_blocs_par_lignes)) {
            cout << "invalid value, (nb_nodes/b * nb_nodes/b) != (nb_blocs_par_lignes * nb_blocs_par_lignes)" << endl; 
            exit(1);
        }
    }
    
    const int K = atoi(argv[2]);


    int *bloc = new int[b*b]();
    
    MPI_Scatter(mat_preparee, b*b, MPI_INT, bloc, b*b, MPI_INT, root, MPI_COMM_WORLD);

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

    int *displs, *sendcount;

    int nb_lignes = (nb_nodes%nprocs >= pid) ? ((nb_nodes/nprocs)+1) : (nb_nodes/nprocs);
    int recvcount = nb_lignes * nb_nodes; 
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
    if (pid == root) {
        cout << "local chosen candidates (reduced)" << endl;
        affichage(reduced_candidates, 1,nb_nodes, 3, INF);
    }

    // TODO NOUVEAUX CODE NATHAN 2025.11.25

    int *medoids = new int[K];

    if (pid == root) {
        //                    max amount of votes a medoid can receive ↓ (if it shoot to high, can add a max() pass before).
        medoids = get_k_best_elt(reduced_candidates, nb_nodes, K, nprocs-1);
    }

    MPI_Bcast(medoids, K, MPI_INT, root, MPI_COMM_WORLD);

    int start_cost = calculate_cost_fragment(medoids, K, mat_distances_fragment, nb_nodes, nb_lignes);
    int global_cost = 0;

    MPI_Reduce(&start_cost, &global_cost, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);

    if (pid == root) {
        // TODO:DELETEME:DEBUG
        cout << "starting cost : " << global_cost << endl ;
        affichage(medoids, 1, K, 3, INF);
    }

    int *permutations = new int[K*K];
    int *permutation_cost = new int[K];
    int *temp = new int[K];


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

    delete[] local_chosen_candidates;
    
    delete[] medoids;
    delete[] permutations;
    delete[] permutation_cost;
    delete[] temp;
    delete[] mat_distances_fragment;
    
    if (pid == root) {
        delete[] mat_adjacence;
        delete[] mat_preparee;
        delete[] mat_gathered;
        delete[] mat_distances;
        delete[] reduced_candidates;
    }    
    
    MPI_Finalize();
    return 0;
}