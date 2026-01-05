#include "needleman_wunsch.hpp"
#include <iostream>

int nw_distance_seq(const char *A, int lenA, const char *B, int lenB) {
	// code inspired by pseudo code presented on 
	// https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm#Advanced_presentation_of_algorithm
	
	int* F = new int [lenA * lenB];

	for (int i = 0; i < lenA; ++i) {
		F[i*lenB] = GAP_PENALTY * i;
	}
	for (int j = 0; j < lenB; ++j) {
		F[j] = GAP_PENALTY * j;
	}
	for (int i = 1; i < lenA; ++i) {
		for (int j = 1; j < lenB; ++j) {
			int _match = F[(i-1) * lenB +  (j-1)] + (A[i] == B[j] ? MATCH_REWARD : MISMATCH_PENALTY);
			int _delete = F[(i-1)*lenB + j] + GAP_PENALTY;
			int _insert = F[i*lenB + (j-1)] + GAP_PENALTY;
			F[i*lenB +j] = max(max(_match, _insert), _delete);
		}
	}
	const int score = F[(lenA * lenB)-1];
	delete[] F;
	return score;
}

int nw_distance_omp(const char *A, int lenA, const char *B, int lenB) {
    return 0;
}
