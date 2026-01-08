#include "Utils.hpp"

void affichage (int* tab, int n, int m, int format, int inf) {
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++)
            // Modified code : affichage du caractère '.' à la place de la valeur inf pour une meilleur lisibilité.
            if (tab[i*m+j] == inf) {
                cout << setw(format) << "." << " ";
            } else {
                cout << setw(format) << tab[i*m+j] << " ";
            }
        cout << endl;
    }
}