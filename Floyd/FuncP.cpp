#include "FuncP.hpp"
#include <iostream>

int *prepareForScatter(int nb_nodes,int *mat_adjacence,int nprocs, int inf)
{   
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
                    // Le bloc ne contient pas la diagonal global, pas de cas particulier traiter.
                    if(!bloc_diagonal){
                        value = inf;
                    }else{
                        // Si x == y alors on est sur un élément de la diagonal global et on garde 0.
                        if(x==y){
                            value = 0;
                        }else{
                            value = inf;
                        }
                    }   
                }
            resultat[indice] = value;
            ++indice;
            }
        }
    }
    return resultat;
}