#pragma once

#define PENALTY -2
#define GAP_PENALTY -1
#define DONE 0
#define UP 1
#define LEFT 2
#define DIAG 3
#define NO_OF_CHARS 256

void matrix_set(int* matrix, int rowSize, int i, int j, int value)
{
	matrix[i * rowSize + j] = value;
}

int matrix_get(int* matrix, int rowSize, int i, int j)
{
	return matrix[i * rowSize + j];
}

double similarityScore(char a, char b)
{
	double result;
	if (a == b)
	{
		result = 2;
	}
	else
	{
		result = PENALTY;
	}
	return result;
}