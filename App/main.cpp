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

    if (argc != 2) {
        cout << "Usage : ./main fichier.dot (graphe au format dot)";
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

        mat_preparee = new int[nb_nodes*nb_nodes]();
        
        mat_preparee = prepareForScatter(nb_nodes, mat_adjacence);
    
        cout << "matrice préparée" << endl;
        affichage(mat_preparee,nb_nodes,nb_nodes,2,INF);
        cout << endl;
    }

    MPI_Bcast(&nb_nodes, 1, MPI_INT, root, MPI_COMM_WORLD);

    int nb_blocs_par_lignes = int(sqrt(double(nprocs)));
    int b = nb_nodes / nb_blocs_par_lignes;

    int *bloc = new int[b*b]();
    
    // if (pid == root) {
    //     int* Dk = MatDistance(nb_nodes, mat_adjacence);
    
    //     cout << "La matrice de distances" << endl;
    //     affichage(Dk,nb_nodes,nb_nodes,3, INF);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);

    
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

    MPI_Barrier(MPI_COMM_WORLD);

    int *mat_gathered, *mat_distances;

    if (pid == root) {
       mat_gathered = new int[nb_nodes*nb_nodes];
       mat_distances = new int[nb_nodes*nb_nodes];
    }

    MPI_Gather(bloc, b*b, MPI_INT, mat_gathered, b*b, MPI_INT, root, MPI_COMM_WORLD);
    
    if (pid == root) { 
        mat_distances = repareAfterGather(nb_nodes, mat_gathered);
        cout << "La matrice de distances" << endl;
        affichage(mat_distances,nb_nodes,nb_nodes,3, INF);
    }
    
    if (pid == root) {
        delete[] mat_adjacence;
        delete[] mat_preparee;
        delete[] mat_gathered;
        delete[] mat_distances;
    }
    delete[] bloc;
    
    MPI_Finalize();
    return 0;
}