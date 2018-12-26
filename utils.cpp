#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include "utils.h"
#include "matrix.h"

void printMatrix(Matrix *matrix) {
	for (int i = 0; i < matrix->size.m; ++i) {
		printf("| ");
		for (int j = 0; j < matrix->size.n; ++j) {
			printf(" %.*f", matrix->matrix[i][j] < 0 ? PRECISION - 1 : PRECISION, matrix->matrix[i][j]);
		}
		printf(" |\n");
	}
	printf("\n");
}

bool cmpSize(Matrix *one, Matrix *two) {
	return one->size.m == two->size.m && one->size.n == two->size.n;
}

bool cmpSizeForMult(Matrix *one, Matrix *two) {
	return one->size.n == two->size.m;
}

double norm(Matrix *matrix) {
	double norm = 0;
	for (int i = 0; i < matrix->size.m; ++i)
		for (int j = 0; j < matrix->size.n; ++j)
			norm += matrix->matrix[i][j] * matrix->matrix[i][j];
	return sqrt(norm);
}

void fillBandMatrix(Matrix *matrix, double *values) {
	for (int i = 0; i < matrix->size.m; ++i) {
		for (int j = 0; j < matrix->size.n; ++j) {
			int diff = abs(i - j);
			if (diff < 3)
				matrix->matrix[i][j] = values[diff];
			else
				matrix->matrix[i][j] = 0;
		}
	}
}

void fillBmatrix(Matrix *matrix) {
	for (int i = 0; i < matrix->size.m; ++i) {
		matrix->matrix[i][0] = sin((i + 1) * 5);
	}
}

void fillMatrix(Matrix *matrix, double value) {
	for (int i = 0; i < matrix->size.m; ++i)
		for (int j = 0; j < matrix->size.n; ++j)
			matrix->matrix[i][j] = value;
}

size_t matrixAllocSize(MatrixSize size) {
	return sizeof(double) * size.m * size.n + sizeof(double *) * size.m;
}

void allocMatrix(Matrix *matrix, bool zero) {
	if (zero)
		matrix->matrix = (double **) calloc(matrixAllocSize(matrix->size), 1);
	else
		matrix->matrix = (double **) malloc(matrixAllocSize(matrix->size));
	if (matrix->matrix == nullptr) {
		printf("failed to alloc matrix");
		return;
	}
	for (int i = 0; i < matrix->size.m; ++i) {
		matrix->matrix[i] = (double *) matrix->matrix + matrix->size.m + i * matrix->size.n;
	}
}

Matrix *createMatrix(MatrixSize size, bool zero) {
	Matrix * matrix = (Matrix*) malloc(sizeof(Matrix));
	matrix->size = size;
	allocMatrix(matrix, zero);
	return matrix;
}

Matrix* createDiagMatrix(MatrixSize size) {
	Matrix *diag = createMatrix(size, true);
	for (int i = 0; i < size.m; ++i)
		diag->matrix[i][i] = 1;
	return diag;
}

void destroyMatrix(Matrix *matrix) {
	free(matrix->matrix);
	free(matrix);
}

Matrix* copyMatrix(Matrix *matrix) {
	Matrix *copy = createMatrix(matrix->size);
	for (int i = 0; i < matrix->size.m; ++i)
		for (int j = 0; j < matrix->size.n; ++j)
			copy->matrix[i][j] = matrix->matrix[i][j];
	return copy;
}
