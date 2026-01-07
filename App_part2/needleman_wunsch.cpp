#include "needleman_wunsch.hpp"

int nw_distance_seq(const char *A, int lenA, const char *B, int lenB)
{
	// code inspired by pseudo code presented on
	// https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm#Advanced_presentation_of_algorithm

	int *F = new int[lenA * lenB];
	bool *H = new bool[lenA * lenB]; // holes

	F[0] = 0;
	H[0] = false;
	for (int i = 1; i < lenA; ++i)
	{
		F[i * lenB] = GAP_PENALTY - i;
		H[i * lenB] = true;
	}
	for (int j = 1; j < lenB; ++j)
	{
		F[j] = GAP_PENALTY - j;
		H[j] = true;
	}
	for (int i = 1; i < lenA; ++i)
	{
		for (int j = 1; j < lenB; ++j)
		{
			_choose_and_write(F, H, i, j, A[i], B[j], lenA, lenB);
		}
	}
	const int score = F[(lenA * lenB) - 1]+1;
	delete[] F;
	return score;
}

int nw_distance_omp(const char *A, int lenA, const char *B, int lenB)
{
	int *F = new int[lenA * lenB];
	bool *H = new bool[lenA * lenB];

	F[0] = 0;

#pragma omp parallel for
	for (int i = 1; i < lenA; ++i)
	{
		F[i * lenB] = GAP_PENALTY - i;
	}

#pragma omp parallel for
	for (int j = 1; j < lenB; ++j)
	{
		F[j] = GAP_PENALTY - j;
	}

#pragma omp parallel
	{
#pragma omp single nowait
		{
			for (int i = 1; i < lenA; ++i)
			{
				for (int j = 1; j < lenB; ++j)
				{
#pragma omp task untied depend(in : F[(i - 1) * lenB + j - 1], F[(i - 1) * lenB + j], F[i * lenB + (j - 1)]) depend(out : F[i * lenB + j])
					{
						_choose_and_write(F, H, i, j, A[i], B[j], lenA, lenB);
					}
				}
			}
		}
	}
	const int score = F[(lenA * lenB) - 1]+1;
	delete[] F;
	return score;
}

void _choose_and_write(int *F, bool *H, int i, int j, char a, char b, int lenA, int lenB)
{
	int _match = F[(i - 1) * lenB + (j - 1)] + (a == b ? MATCH_REWARD : MISMATCH_PENALTY);
	int _delete = F[(i - 1) * lenB + j] + ((H[(i - 1) * lenB + j]) ? EXTEND_PENALTY : GAP_PENALTY);
	int _insert = F[i * lenB + (j - 1)] + ((H[i * lenB + (j - 1)]) ? EXTEND_PENALTY : GAP_PENALTY);
	int choice = max(max(_match, _insert), _delete);
	F[i * lenB + j] = choice;
	H[i * lenB + j] = choice != _match;
}
