/**********************************************
File		:	algebra.cpp
Author		:	Mingcheng Chen
Last Update	:	January 14th, 2014
***********************************************/

#include "algebra.h"
#include "utility.h"

lcs::Matrix lcs::MatrixMatrixMultiplication(const Matrix &a, const Matrix &b) {
    if (a.numOfCols != b.numOfRows)
        lcs::Error("Matrix sizes do not match.");

    lcs::Matrix c(a.numOfRows, b.numOfCols);
    for (int i = 0; i < a.numOfRows; i++)
        for (int j = 0; j < b.numOfCols; j++)
            for (int k = 0; k < a.numOfCols; k++)
                MAT(c, i, j) += MAT(a, i, k) * MAT(b, k, j);

    return c;
}

lcs::Matrix lcs::MatrixTranspose(const Matrix &a) {
    lcs::Matrix b(a.numOfCols, a.numOfRows);
    for (int i = 0; i < a.numOfRows; i++)
        for (int j = 0; j < a.numOfCols; j++)
            MAT(b, j, i) = MAT(a, i, j);
    return b;
}
