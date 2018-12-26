#include <cstdlib>
#include <cstdio>
#include "matrix.h"
#include "utils.h"

#define MATRIX_DIM 996

int main() {
	Matrix *aMatrix = createMatrix({MATRIX_DIM, MATRIX_DIM});
	double values[] = {11.0, -1.0, -1.0};
	fillBandMatrix(aMatrix, values);
	Matrix *bMatrix = createMatrix({MATRIX_DIM, 1});
	fillBmatrix(bMatrix);
	double precision = 10e-9;

	printf("ZADANIE B\n");
	solveIterative(aMatrix, bMatrix, precision, JACOBI);
	solveIterative(aMatrix, bMatrix, precision, GAUSS_SEIDEL);

	printf("\nZADANIE C\n");
	values[0] = 3;
	fillBandMatrix(aMatrix, values);
	solveIterative(aMatrix, bMatrix, precision, JACOBI);
	solveIterative(aMatrix, bMatrix, precision, GAUSS_SEIDEL);


	printf("\nZADANIE D\n");
	solveLU(aMatrix, bMatrix);

	printf("\nZADANIE E\n");
	int N[] = {100, 500, 1000, 2000, 3000, 4000, 5000};
	int nSize = 7;
	printf("N Jacobi Gauss-Seidel LU\n");
	double time;
	values[0] = 11;
	for (int i = 0; i < nSize; ++i) {
		destroyMatrix(aMatrix);
		destroyMatrix(bMatrix);
		aMatrix = createMatrix({N[i], N[i]});
		fillBandMatrix(aMatrix, values);
		bMatrix = createMatrix({N[i], 1});
		fillBmatrix(bMatrix);
		printf("%d ", N[i]);
		solveIterative(aMatrix, bMatrix, precision, JACOBI, false, false, &time);
		printf("%fs ", time);
		solveIterative(aMatrix, bMatrix, precision, GAUSS_SEIDEL, false, false, &time);
		printf("%fs ", time);
		solveLU(aMatrix, bMatrix, false, false, &time);
		printf("%fs\n", time);
	}
	destroyMatrix(aMatrix);
	destroyMatrix(bMatrix);
	return 0;
}
