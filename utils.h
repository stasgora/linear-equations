#ifndef II_UTILS_H
#define II_UTILS_H

#include "matrix.h"

#define PRECISION 6

void printMatrix(Matrix *matrix);

bool cmpSize(Matrix *one, Matrix *two);

bool cmpSizeForMult(Matrix *one, Matrix *two);

double norm(Matrix *matrix);

void fillBandMatrix(Matrix *matrix, double *values);

void fillBmatrix(Matrix *matrix);

void fillMatrix(Matrix* matrix, double value);

Matrix* createDiagMatrix(MatrixSize size);

double** allocMatrix(MatrixSize size, bool zero = false);

Matrix* createMatrix(MatrixSize size, bool zero = false);

void destroyMatrix(Matrix *matrix);

Matrix* copyMatrix(Matrix* matrix);

#endif //II_UTILS_H
