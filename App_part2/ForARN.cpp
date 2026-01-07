#include "ForARN.hpp"
#include <mpi.h>

char *readArnFromFile(char *f, int nb_nodes) {
	char* out = new char[(nb_nodes * TAILLESEQ)];
	int cursor = 0;
	string currentline;
	ifstream File(f);

	while(getline(File, currentline)) {
		// exemple
		// >2
		// TGCCCTTCCTCTTCCTCGTAAGGAAGTAGACTACATTACGGGAGGCCTTAAGTAGCCTTAAATTCAGGCA
		// AGCATCGTTCTTTGAGTTGGTTGCACCCAA

		// Si la ligne commence par '>', on l'ignore
		if (currentline[0] != '>') {
			// sinon
			// nombre de charactère de la ligne courante
			int nb_char = currentline.length();
			// copy dans notre tableau la ligne courante
			memcpy(&out[cursor], currentline.c_str() , nb_char * sizeof(char));
			// on décale la tête de lecture de notre tableau
			cursor += nb_char;
		}
	}

	return out ;
}

void buildMatrixFragment(int *mat_adjacence, int nb_nodes, char *sequence_fragment, int fragment_offset, int nb_lignes, char *recv_buffer, int recv_offset, int recv_nb_lignes) {
	// mat_adjacence est seulement modifiée ici, pas besoin de malloc de la mémoire
	int distance;
	// pour chaque séquence d'ARN dans sequence_fragment
	for (int i = 0; i < nb_lignes; ++i) {
		// pour chaque séquence d'ARN dans recv_buffer
		for (int j = 0; j < recv_nb_lignes; ++j) {
			distance = nw_distance_omp(&sequence_fragment[i*TAILLESEQ], TAILLESEQ, &recv_buffer[j*TAILLESEQ], TAILLESEQ);
			if (distance > 0) {
				mat_adjacence[(fragment_offset+i) *nb_nodes + (recv_offset+j)] = TAILLESEQ - distance;
			}
			// sinon on laisse à 0
		}
	}
}
