#ifndef NEEDLEMAN
#define NEEDELMAN

#include <algorithm> // needed for max()
#include <omp.h>

using namespace std;

#define MATCH_REWARD 1

#define MISMATCH_PENALTY -1

#define GAP_PENALTY -3

int nw_distance_seq(const char *A, int lenA, const char *B, int lenB);

int nw_distance_omp(const char *A, int lenA, const char *B, int lenB);

void _choose_and_write(int* F, int i, int j, char a, char b, int lenA, int lenB);

#endif