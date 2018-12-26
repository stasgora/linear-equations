#ifndef II_OPERATIONS_H
#define II_OPERATIONS_H

#define JACOBI 0
#define GAUSS_SEIDEL 1
#define ITER_LIMIT 10000

typedef struct {
	int m, n;
} MatrixSize;

typedef struct {
	MatrixSize size;
	double **matrix;
} Matrix;

void divideIntoLUD(Matrix *A, Matrix **L, Matrix **U, Matrix **D);

Matrix *solveIterative(Matrix* A, Matrix *b, double precision, int method,
					   bool returnX = false, bool printStats = true, double *time = nullptr);

Matrix *add(Matrix *one, Matrix *two, double mult = 1, bool intoOne = false);

Matrix *mult(Matrix *one, Matrix *two, Matrix *res = nullptr);

void LUfactor(Matrix *A, Matrix **L, Matrix **U);

Matrix *solveLU(Matrix *A, Matrix *b, bool returnX = false, bool printStats = true, double *time = nullptr);

Matrix *substitute(Matrix *A, Matrix *b, bool front = true, Matrix *res = nullptr);

#endif //II_OPERATIONS_H
