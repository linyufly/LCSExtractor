/**********************************************
File		:	algebra.cpp
Author		:	Mingcheng Chen
Last Update	:	January 14th, 2014
***********************************************/

#include "algebra.h"
#include "utility.h"

#include <vtkSmartPointer.h>
#include <vtkPCAStatistics.h>
#include <vtkDoubleArray.h>
#include <vtkTable.h>
#include <vtkIndent.h>
#include <vtkMath.h>

#include <iostream>

lcs::Matrix lcs::MatrixMatrixMultiplication(const lcs::Matrix &a, const lcs::Matrix &b) {
    if (a.numOfCols != b.numOfRows)
        lcs::Error("Matrix sizes do not match.");

    lcs::Matrix c(a.numOfRows, b.numOfCols);
    for (int i = 0; i < a.numOfRows; i++)
        for (int j = 0; j < b.numOfCols; j++)
            for (int k = 0; k < a.numOfCols; k++)
                MAT(c, i, j) += MAT(a, i, k) * MAT(b, k, j);

    return c;
}

lcs::Matrix lcs::MatrixTranspose(const lcs::Matrix &a) {
    lcs::Matrix b(a.numOfCols, a.numOfRows);
    for (int i = 0; i < a.numOfRows; i++)
        for (int j = 0; j < a.numOfCols; j++)
            MAT(b, j, i) = MAT(a, i, j);
    return b;
}

lcs::Matrix lcs::PCA(const lcs::Matrix &a) {
    Matrix C = lcs::MatrixMatrixMultiplication(lcs::MatrixTranspose(a), a);
    double **CForVTK = C.MatrixForVTK();
    double eigenvalues[3];
    double **eigenvectors = lcs::Matrix::CreateVTKMatrix(3, 3);
    vtkMath::Jacobi(CForVTK, eigenvalues, eigenvectors);
    delete [] CForVTK;
    lcs::Matrix result(eigenvectors[0], 3, 3);
    lcs::Matrix::DisposeVTKMatrix(eigenvectors);
    return result;
/*
    static char str[100];

    vtkSmartPointer<vtkDoubleArray> *data = new vtkSmartPointer<vtkDoubleArray> [a.numOfCols];
    for (int i = 0; i < a.numOfCols; i++) {
        data[i] = vtkSmartPointer<vtkDoubleArray>::New();
        sprintf(str, "C%d", i);
        data[i]->SetName(str);
        data[i]->SetNumberOfComponents(1);
        for (int j = 0; j < a.numOfRows; j++)
            data[i]->InsertNextValue(MAT(a, j, i));
    }

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
    for (int i = 0; i < a.numOfCols; i++)
        table->AddColumn(data[i]);

    vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();
    pcaStatistics->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, table);

    for (int i = 0; i < a.numOfCols; i++) {
        sprintf(str, "C%d", i);
        pcaStatistics->SetColumnStatus(str, 1);
    }
    pcaStatistics->RequestSelectedColumns();
    pcaStatistics->SetDeriveOption(true);
    pcaStatistics->Update();

    /// DEBUG ///
    //vtkSmartPointer<vtkDoubleArray> eigenvalues = vtkSmartPointer<vtkDoubleArray>::New();
    //pcaStatistics->GetEigenvalues(eigenvalues);

    //for (vtkIdType i = 0; i < eigenvalues->GetNumberOfTuples(); i++)
    //    printf("Eigenvalue %d = %lf\n", (int)i, eigenvalues->GetValue(i));

    vtkSmartPointer<vtkDoubleArray> eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
    pcaStatistics->GetEigenvectors(eigenvectors);

    lcs::Matrix result(a.numOfRows, a.numOfCols);

    /// DEBUG ///
    //vtkIndent indent(0);
    double *evec = new double [eigenvectors->GetNumberOfComponents()];
    for (vtkIdType i = 0; i < eigenvectors->GetNumberOfTuples(); i++) {
        //printf("Eigenvector %d: ", (int)i);
       eigenvectors->GetTuple(i, evec);
        for (vtkIdType j = 0; j < eigenvectors->GetNumberOfComponents(); j++) {
            //printf("%lf ", evec[j]);
            MAT(result, i, j) = evec[j];
        }

        /// DEBUG ///
        //printf("|");

        //vtkSmartPointer<vtkDoubleArray> eigenvectorSingle = vtkSmartPointer<vtkDoubleArray>::New();
        //pcaStatistics->GetEigenvector(i, eigenvectorSingle);
        //eigenvectorSingle->GetTuple(0, evec);
        //for (vtkIdType j = 0; j < eigenvectorSingle->GetNumberOfComponents(); j++)
        //    printf(" %lf", evec[j]);
        //eigenvectorSingle->PrintSelf(std::cout, indent);
        //printf("\n");
    }
    delete [] evec;

    return result;
*/
}
