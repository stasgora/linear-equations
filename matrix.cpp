#include <cstdio>
#include <cfloat>
#include <chrono>
#include <ctime>
#include <cmath>
#include "matrix.h"
#include "utils.h"

void divideIntoLUD(Matrix *A, Matrix **L, Matrix **U, Matrix **D) {
	*L = createMatrix(A->size, true);
	*U = createMatrix(A->size, true);
	*D = createMatrix(A->size, true);
	for (int i = 0; i < A->size.m; ++i) {
		for (int j = 0; j < A->size.n; ++j) {
			if(i < j)
				(*U)->matrix[i][j] = -A->matrix[i][j];
			else if(i > j)
				(*L)->matrix[i][j] = -A->matrix[i][j];
			else
				(*D)->matrix[i][j] = A->matrix[i][j];
		}
	}
}

Matrix *solveIterative(Matrix *A, Matrix *b, double precision, int method, bool returnX, bool printStats, double *time) {
	clock_t start = clock();

	Matrix *L, *U, *D, *x;
	divideIntoLUD(A, &L, &U, &D);

	x = createMatrix(b->size);
	fillMatrix(x, 1.0 / A->size.m);

	Matrix *constSum, *constSubst;
	if(method == JACOBI) {
		constSum = add(L, U);
		constSubst = substitute(D, b, false);
	} else {
		constSum = add(D, L, -1);
		constSubst = substitute(constSum, b);
	}
	Matrix *temp = createMatrix(b->size);
	double residuum = DBL_MAX;
	int iterations = 0;

	while (residuum > precision) {
		if (iterations > ITER_LIMIT)
			break;
		if(method == JACOBI) {
			mult(constSum, x, temp);
			substitute(D, temp, false, x);
		} else {
			mult(U, x, temp);
			substitute(constSum, temp, true, x);
		}
		add(x, constSubst, 1, true);

		mult(A, x, temp);
		add(temp, b, -1, true);
		residuum = norm(temp);
		iterations++;
	}
	destroyMatrix(L);
	destroyMatrix(U);
	destroyMatrix(D);
	destroyMatrix(constSum);
	destroyMatrix(constSubst);
	destroyMatrix(temp);

	double runtime = ((double) clock() - start) / CLOCKS_PER_SEC;
	if(time != nullptr)
		*time = runtime;
	if(printStats) {
		printf(method == JACOBI ? " -Metoda Jacobiego\n" : " -Metoda Gaussa-Seidla\n");

		if (iterations <= ITER_LIMIT) {
			printf("\tCzas działania: %fs\n", runtime);
			printf("\tIlość iteracji: %d\n", iterations);
			//printf("\tNorma residuum: %e\n", residuum);

			bool NaN = false, inf = false;
			for (int i = 0; i < x->size.m; ++i)
				if(!std::isfinite(x->matrix[i][0])) {
					if (!inf && std::isinf(x->matrix[i][0])) {//std::isnan(x->matrix[i][0])) {
						printf("\tWykryto zbieżność do: %f - nieskończoności\n", x->matrix[i][0]);
						inf = true;
					} else if(!NaN) {
						printf("\tWykryto wartość %f - nie-numer\n", x->matrix[i][0]);
						NaN = true;
					}
					if (NaN && inf)
						break;
				}
		} else
			printf("\tBrak zbieżności\n");
	}
	if (returnX)
		return x;
	destroyMatrix(x);
	return nullptr;
}

Matrix* add(Matrix *one, Matrix *two, double mult, bool intoOne) {
	if (!cmpSize(one, two)) {
		printf("matrix sizes don't match\n");
		return nullptr;
	}
	Matrix *res = nullptr;
	if (!intoOne)
		res = createMatrix(one->size);
	for (int i = 0; i < one->size.m; ++i)
		for (int j = 0; j < one->size.n; ++j)
			if (!intoOne)
				res->matrix[i][j] = one->matrix[i][j] + mult * two->matrix[i][j];
			else
				one->matrix[i][j] = one->matrix[i][j] + mult * two->matrix[i][j];
	if (!intoOne)
		return res;
	else
		return one;
}

Matrix *mult(Matrix *one, Matrix *two, Matrix *res) {
	if (!cmpSizeForMult(one, two) || (res != nullptr && (one->size.m != res->size.m || two->size.n != res->size.n))) {
		printf("matrix sizes don't match\n");
		return nullptr;
	}
	if (res == nullptr)
		res = createMatrix({one->size.m, two->size.n});
	for (int i = 0; i < one->size.m; ++i)
		for (int j = 0; j < two->size.n; ++j) {
			double sum = 0;
			for (int k = 0; k < one->size.n; ++k)
				sum += one->matrix[i][k] * two->matrix[k][j];
			res->matrix[i][j] = sum;
		}
	return res;
}

Matrix *solveLU(Matrix *A, Matrix *b, bool returnX, bool printStats, double *time) {
	clock_t start = clock();

	Matrix *L, *U;
	LUfactor(A, &L, &U);
	Matrix *y = substitute(L, b);
	Matrix *x = substitute(U, y, false);
	destroyMatrix(L);
	destroyMatrix(U);

	mult(A, x, y);
	add(y, b, -1, true);
	double residuum = norm(y);
	destroyMatrix(y);

	double runtime = ((double) clock() - start) / CLOCKS_PER_SEC;
	if(time != nullptr)
		*time = runtime;
	if(printStats) {
		printf(" -Faktoryzacja LU\n");
		printf("\tCzas działania: %fs\n", runtime);
		printf("\tNorma residuum: %e\n", residuum);
	}
	if (returnX)
		return x;
	destroyMatrix(x);
	return nullptr;
}

void LUfactor(Matrix *A, Matrix **L, Matrix **U) {
	*U = copyMatrix(A);
	*L = createDiagMatrix(A->size);
	for (int k = 0; k < A->size.m - 1; ++k) {
		for (int j = k + 1; j < A->size.m; ++j) {
			(*L)->matrix[j][k] = (*U)->matrix[j][k] / (*U)->matrix[k][k];
			for (int i = k; i < A->size.m; ++i)
				(*U)->matrix[j][i] -= (*L)->matrix[j][k] * (*U)->matrix[k][i];
		}
	}
}

Matrix *substitute(Matrix *A, Matrix *b, bool front, Matrix *res) {
	int start = front ? 0 : A->size.m - 1, end = front ? A->size.m : 0;
	if (res == nullptr)
		res = createMatrix(b->size);
	for (int i = start; front ? i < end : i >= end; i += front ? 1 : -1) {
		double sum = b->matrix[i][0];
		int start2 = front ? 0 : i + 1, end2 = front ? i : A->size.m;
		for (int j = start2; j < end2; j++)
			sum -= res->matrix[j][0] * A->matrix[i][j];
		res->matrix[i][0] = sum / A->matrix[i][i];
	}
	return res;
}
