#include "FuncP.hpp"

int *prepareForScatter(int nb_nodes,int *mat_adjacence,int nprocs)
{   
    int nb_blocs_par_lignes = int(sqrt(double(nprocs)));
    int b = nb_nodes / nb_blocs_par_lignes;

    // matrice envoyée
    int* resultat = new int[nb_nodes*nb_nodes]();

    // position de la "tête d'écriture" dans résultat
    int indice = 0;

    // On parcours les blocs, (chaque processeur reçevra un blocs, nb_blocs = nprocs)
    for (int bloc_id = 0; bloc_id < nprocs; ++bloc_id) {
        int numero_colonne_bloc = bloc_id % nb_blocs_par_lignes;
        int numero_ligne_bloc = bloc_id / nb_blocs_par_lignes;

        // position du premier element du bloc dans la matrice d'adjacence
        int position_debut_bloc = (numero_ligne_bloc * (b*b) * nb_blocs_par_lignes) + (numero_colonne_bloc * b);

        // Pour chaque element du bloc
        for (int x = 0; x < b; ++x) {
            for (int y = 0; y < b; ++y) {
                // on ecrit l'element dans le résultat
                // + (x * nb_nodes) permet de se décaler d'une ligne vers le bas dans la matrice d'adjacence.
                resultat[indice] = mat_adjacence[position_debut_bloc + (x * nb_nodes) + y];
                ++indice;
            }
        }
    }
    return resultat;
}