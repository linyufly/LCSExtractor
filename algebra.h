/**********************************************
File		:	algebra.h
Author		:	Mingcheng Chen
Last Update	:	January 14th, 2014
***********************************************/

#ifndef __ALGEBRA_H
#define __ALGEBRA_H

#include <cstdio>
#include <cstring>

#define MAT(mat, i, j) (mat).arr[(i) * (mat).numOfCols + (j)]

namespace lcs {

class Matrix {
private:
    double *arr;
    int numOfRows, numOfCols;

public:
    Matrix() {
        arr = NULL;
        numOfRows = numOfCols = 0;
    }

    ~Matrix() {
        if (arr) delete [] arr;
    }

    Matrix(int numOfRows, int numOfCols) {
        this->arr = new double [numOfRows * numOfCols];
        this->numOfRows = numOfRows;
        this->numOfCols = numOfCols;
        memset(this->arr, 0, sizeof(double) * numOfRows * numOfCols);
    }

    Matrix(double *arr, int numOfRows, int numOfCols) {
        this->arr = new double [numOfRows * numOfCols];
        this->numOfRows = numOfRows;
        this->numOfCols = numOfCols;
        memcpy(this->arr, arr, sizeof(double) * numOfRows * numOfCols);
    }

    Matrix &operator = (const Matrix &a) {
        if (this == &a) return *this;

        this->numOfRows = a.numOfRows;
        this->numOfCols = a.numOfCols;
        if (this->arr) delete this->arr;
        this->arr = new double [a.numOfRows * a.numOfCols];
        memcpy(this->arr, a.arr, sizeof(double) * a.numOfRows * a.numOfCols);

        return *this;
    }

    double **MatrixForVTK() const {
        double **mat = new double * [this->numOfRows];
        for (int i = 0; i < this->numOfRows; i++)
            mat[i] = this->arr + i * this->numOfCols;
        return mat;
    }

    int GetNumOfRows() const {
        return numOfRows;
    }

    int GetNumOfCols() const {
        return numOfCols;
    }

    void Output() const {
        printf("#BEGIN\n");
        for (int i = 0; i < numOfRows; i++) {
            for (int j = 0; j < numOfCols; j++)
                printf(" %lf", arr[i * numOfCols + j]);
            printf("\n");
        }
        printf("#END\n");
    }

    double &Element(int row, int col) const {
        return MAT(*this, row, col);
    }

    static double **CreateVTKMatrix(int numOfRows, int numOfCols) {
        double *arr = new double [numOfRows * numOfCols];
        double **mat = new double * [numOfRows];
        for (int i = 0; i < numOfRows; i++)
            mat[i] = arr + i * numOfCols;
        return mat;
    }

    static void DisposeVTKMatrix(double **mat) {
        delete [] mat[0];
        delete [] mat;
    }

    friend Matrix MatrixMatrixMultiplication(const Matrix &a, const Matrix &b);
    friend Matrix MatrixTranspose(const Matrix &a);
    friend Matrix PCA(const Matrix &a);
};

Matrix MatrixMatrixMultiplication(const Matrix &a, const Matrix &b);
Matrix MatrixTranspose(const Matrix &a);
Matrix PCA(const Matrix &a);

}

#endif
