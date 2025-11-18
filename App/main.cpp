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

        cout << "matrice d'adjacence" << endl;
        affichage(mat_adjacence,nb_nodes,nb_nodes,2, INF);
        cout << endl;
        
        mat_preparee = prepareForScatter(nb_nodes, mat_adjacence);
        
        cout << "matrice préparée" << endl;
        affichage(mat_preparee,nb_nodes,nb_nodes,2,INF);
        cout << endl;
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

    // if we need to shuffle, do it here.

    int *displs, *sendcount;

    int recvcount = (nb_nodes%nprocs >= pid) ? (((nb_nodes/nprocs)+1)*nb_nodes) : ((nb_nodes/nprocs)*nb_nodes); 
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
    // On a plus besoin du tableau des candidats locaux
    delete[] local_chosen_candidates;

    std::vector<std::vector<int>>* candidates_per_value;
    int* candidates_globaux;
    int* cout_globaux_reduced;
    if(pid != root){
        //Buffer de reception lors du broadcast des candidats.
        candidates_globaux = new int[nb_nodes];
    }
    if (pid == root) {
        cout_globaux_reduced = new int[nb_nodes*K]; // *K ?
        candidates_per_value = new std::vector<std::vector<int>>;
        // nprocs+1 car certains processus peuvent "voter" pour un noeud n'étant pas dans leur fragment, et donc les valeurs peuvent aller de 0 à nprocs inclus.
        (*candidates_per_value).resize(nprocs+1); // retirer le -1 une fois l'init fix
        affichage(reduced_candidates, 1, nb_nodes, 2, INF);
        candidates_globaux = process_candidates(candidates_per_value,reduced_candidates,nb_nodes,K);
        affichage(candidates_globaux, 1, K, 2, INF);
        delete candidates_per_value;
    }

    // On écrira dans ce tableau le cout local pour un changement de médoide pour un autre noeud. Et ce, pour l'ensemble des médoides (donc K tableau de cout).

    int* cout_locaux = new int[nb_nodes*K];
    int changement = 1;
    while(changement){
        changement = 0;
        MPI_Bcast(candidates_globaux,nb_nodes,MPI_INT,root,MPI_COMM_WORLD);
        calcul_cout_swap(candidates_globaux, cout_locaux,mat_distances_fragment, K, nb_nodes, recvcount/nb_nodes);
        MPI_Reduce(cout_locaux,cout_globaux_reduced,nb_nodes*K,MPI_INT,MPI_SUM,root,MPI_COMM_WORLD);
        if(pid == root){
            changement = choix_nouveaux_candidats(nb_nodes,cout_globaux_reduced, K,candidates_globaux);
        }
        MPI_Bcast(&changement,1,MPI_INT,root,MPI_COMM_WORLD);
    }

    if (pid == root) {
        affichage(candidates_globaux, 1, nb_nodes, 2, INF);
    }    



    delete[] cout_locaux;
    delete[] mat_distances_fragment;
    
    if (pid == root) {
        delete[] cout_globaux_reduced;
        delete[] candidates_globaux;
        delete[] mat_adjacence;
        delete[] mat_preparee;
        delete[] mat_gathered;
        delete[] mat_distances;
        delete[] reduced_candidates;
    }
    
    MPI_Finalize();
    return 0;
}