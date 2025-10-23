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
        affichage(mat_adjacence,nb_nodes,nb_nodes,2);
        cout << endl;

        mat_preparee = new int[nb_nodes*nb_nodes]();
        
        mat_preparee = prepareForScatter(nb_nodes, mat_adjacence, nprocs);
    
        cout << "matrice préparée" << endl;
        affichage(mat_preparee,nb_nodes,nb_nodes,2);
        cout << endl;
    }

    MPI_Bcast(&nb_nodes, 1, MPI_INT, root, MPI_COMM_WORLD);

    int nb_blocs_par_lignes = int(sqrt(double(nprocs)));
    int b = nb_nodes / nb_blocs_par_lignes;

    int *bloc = new int[b*b]();
    
    // int* Dk = MatDistance(nb_nodes, mat_adjacence);

    // cout << "La matrice de distances" << endl;
    // affichage(Dk,nb_nodes,nb_nodes,3);
    
    MPI_Scatter(mat_preparee, b*b, MPI_INT, bloc, b*b, MPI_INT, root, MPI_COMM_WORLD);

    if (pid == 0) {
        cout << "matrice d'envoyé à 0" << endl;
        affichage(bloc,b,b,2);
        cout << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (pid == root) {
        delete[] mat_adjacence;
        delete[] mat_preparee;
    }
    delete[] bloc;

    MPI_Finalize();
    return 0;
}