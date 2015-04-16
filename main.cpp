/**********************************************
File		:	main.cpp
Author		:	Mingcheng Chen
Last Update	:	January 21st, 2014
***********************************************/

#include "lcsGeometry.h"
#include "utility.h"
#include "algebra.h"
#include "marchingCubesTable.h"

#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <list>
#include <algorithm>

lcs::Configuration configure;
lcs::Vector *lastPositions;
double *ftleValues;
lcs::Vector *ftleGradients;
lcs::Vector *N1;

void ReadConfiguration() {
    configure.LoadFile("configuration.conf");
    printf("Finished reading the configuration file.\n");
    printf("\n");
}

void ReadLastPositions() {
    printf("Start reading the data file ...\n");

    int nx = configure.GetNx();
    int ny = configure.GetNy();
    int nz = configure.GetNz();
    int length = (nx + 1) * (ny + 1) * (nz + 1);
    lastPositions = new lcs::Vector [length];

    FILE *fin = fopen(configure.GetDataFile().c_str(), "r");
    if (fin == NULL)
        lcs::Error("The data file \"%s\" cannot be found.", configure.GetDataFile().c_str());

    for (int i = 0; i < length; i++) {
        int a, b, c;
        double x, y, z;
        fscanf(fin, "%d %d %d: %lf %lf %lf", &a, &b, &c, &x, &y, &z);
        lastPositions[lcs::Code(a, b, c, ny, nz)] = lcs::Vector(x, y, z);
    }

    printf("Finished reading the data file.\n");
    printf("\n");
}

#define P(i, j, k) lastPositions[lcs::Code((i), (j), (k), ny, nz)]
#define FTLE(i, j, k) ftleValues[lcs::Code((i), (j), (k), ny, nz)]
#define N1(i, j, k) N1[lcs::Code((i), (j), (k), ny, nz)]

void OutputImageData(vtkImageData* data, int nx, int ny, int nz, const char *fileName) {
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(fileName);
    writer->SetInputData(data);
    writer->Write();
}

void GetRawFTLE() {
    printf("Start to get the raw FTLE ...\n");

    int nx = configure.GetNx();
    int ny = configure.GetNy();
    int nz = configure.GetNz();
    int stencilSize = configure.GetStencilSize();
    int length = (nx + 1) * (ny + 1) * (nz + 1);
    int delta = stencilSize / 2;

    double dx = configure.GetDx();
    double dy = configure.GetDy();
    double dz = configure.GetDz();

    ftleValues = new double [length];
    N1 = new lcs::Vector [length];

    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k <= nz; k++) {
                double dx_x, dy_x, dz_x,
                       dx_y, dy_y, dz_y,
                       dx_z, dy_z, dz_z;
                // dx
                if (!i) {
                    dx_x = (P(i + 1, j, k).GetX() - P(i, j, k).GetX()) / dx;
                    dx_y = (P(i + 1, j, k).GetY() - P(i, j, k).GetY()) / dx;
                    dx_z = (P(i + 1, j, k).GetZ() - P(i, j, k).GetZ()) / dx;
                } else if (i == nx) {
                    dx_x = (P(i, j, k).GetX() - P(i - 1, j, k).GetX()) / dx;
                    dx_y = (P(i, j, k).GetY() - P(i - 1, j, k).GetY()) / dx;
                    dx_z = (P(i, j, k).GetZ() - P(i - 1, j, k).GetZ()) / dx;
                } else {
                    dx_x = (P(i + 1, j, k).GetX() - P(i - 1, j, k).GetX()) / (2 * dx);
                    dx_y = (P(i + 1, j, k).GetY() - P(i - 1, j, k).GetY()) / (2 * dx);
                    dx_z = (P(i + 1, j, k).GetZ() - P(i - 1, j, k).GetZ()) / (2 * dx);
                }
                // dy
                if (!j) {
                    dy_x = (P(i, j + 1, k).GetX() - P(i, j, k).GetX()) / dy;
                    dy_y = (P(i, j + 1, k).GetY() - P(i, j, k).GetY()) / dy;
                    dy_z = (P(i, j + 1, k).GetZ() - P(i, j, k).GetZ()) / dy;
                } else if (j == ny) {
                    dy_x = (P(i, j, k).GetX() - P(i, j - 1, k).GetX()) / dy;
                    dy_y = (P(i, j, k).GetY() - P(i, j - 1, k).GetY()) / dy;
                    dy_z = (P(i, j, k).GetZ() - P(i, j - 1, k).GetZ()) / dy;
                } else {
                    dy_x = (P(i, j + 1, k).GetX() - P(i, j - 1, k).GetX()) / (2 * dy);
                    dy_y = (P(i, j + 1, k).GetY() - P(i, j - 1, k).GetY()) / (2 * dy);
                    dy_z = (P(i, j + 1, k).GetZ() - P(i, j - 1, k).GetZ()) / (2 * dy);
                }
                // dz
                if (!k) {
                    dz_x = (P(i, j, k + 1).GetX() - P(i, j, k).GetX()) / dz;
                    dz_y = (P(i, j, k + 1).GetY() - P(i, j, k).GetY()) / dz;
                    dz_z = (P(i, j, k + 1).GetZ() - P(i, j, k).GetZ()) / dz;
                } else if (k == nz) {
                    dz_x = (P(i, j, k).GetX() - P(i, j, k - 1).GetX()) / dz;
                    dz_y = (P(i, j, k).GetY() - P(i, j, k - 1).GetY()) / dz;
                    dz_z = (P(i, j, k).GetZ() - P(i, j, k - 1).GetZ()) / dz;
                } else {
                    dz_x = (P(i, j, k + 1).GetX() - P(i, j, k - 1).GetX()) / (2 * dz);
                    dz_y = (P(i, j, k + 1).GetY() - P(i, j, k - 1).GetY()) / (2 * dz);
                    dz_z = (P(i, j, k + 1).GetZ() - P(i, j, k - 1).GetZ()) / (2 * dz);
                }

                double arr[] = {dx_x, dy_x, dz_x, dx_y, dy_y, dz_y, dx_z, dy_z, dz_z};

                /// DEBUG ///
                if (i == 10 && j == 20 && k == 30) {
                  printf("f\n");
                  for (int row = 0; row < 3; row++) {
                    for (int col = 0; col < 3; col++) {
                      printf("%lf ", arr[row * 3 + col]);
                    }
                    printf("\n");
                  }
                }

                lcs::Matrix F(arr, 3, 3);
                lcs::Matrix C = lcs::MatrixMatrixMultiplication(lcs::MatrixTranspose(F), F);

                double eigenvalues[3];
                double **eigenvectors = lcs::Matrix::CreateVTKMatrix(3, 3);
                double **CForVTK = C.MatrixForVTK();

                /// DEBUG ///
                if (i == 10 && j == 20 && k == 30) {
                  printf("cc\n");
                  for (int row = 0; row < 3; row++) {
                    for (int col = 0; col < 3; col++) {
                      printf("%lf ", CForVTK[row][col]);
                    }
                    printf("\n");
                  }
                }

                vtkMath::Jacobi(CForVTK, eigenvalues, eigenvectors);
                delete [] CForVTK;

                /// DEBUG ///
                if (i == 10 && j == 20 && k == 30) {
                  printf("eigenvalue = %lf\n", eigenvalues[0]);
                }

                FTLE(i, j, k) = log(sqrt(eigenvalues[0]));
                N1(i, j, k) = lcs::Vector(eigenvectors[0][0], eigenvectors[1][0], eigenvectors[2][0]);

                lcs::Matrix::DisposeVTKMatrix(eigenvectors);
            }

    vtkSmartPointer<vtkImageData> image3D = vtkSmartPointer<vtkImageData>::New();
    image3D->SetExtent(0, nx, 0, ny, 0, nz);
    image3D->AllocateScalars(VTK_DOUBLE, 1);

    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k <= nz; k++)
                image3D->SetScalarComponentFromDouble(i, j, k, 0, FTLE(i, j, k));

    OutputImageData(image3D, nx, ny, nz, "RawFTLEvalues.vtu");

    printf("Done.\n\n");
}

double GaussianGradientCoefficient(double sigma) {
    return -1.0 / lcs::Sqr(sigma) / (sqrt(lcs::Cub(2 * vtkMath::Pi())) * lcs::Cub(sigma)); 
}

double GaussianGradientExponent(double n, double sigma) {
    return exp(-lcs::Sqr(n) / (2 * lcs::Sqr(sigma)));
}

#define _F(i, j, k) F[lcs::Code((i), (j), (k), ny + delta * 2, nz + delta * 2)]
#define dxF(i, j, k) dxF[lcs::Code((i), (j), (k), ny, nz)]
#define dyF(i, j, k) dyF[lcs::Code((i), (j), (k), ny, nz)]
#define dzF(i, j, k) dzF[lcs::Code((i), (j), (k), ny, nz)]
#define S1(i, j, k) stage1[lcs::Code((i), (j), (k), ny + delta * 2, nz + delta * 2)]
#define S2(i, j, k) stage2[lcs::Code((i), (j), (k), ny + delta * 2, nz + delta * 2)]

void GetSmoothedDerivative(double *F, double *dxF, double *dyF, double *dzF, int nx, int ny, int nz, double dx, double dy, double dz, double sigma, int stencilSize) {
    int length = (nx + stencilSize) * (ny + stencilSize) * (nz + stencilSize);
    int delta = stencilSize / 2;
    
    double *stage1 = new double [length];
    double *stage2 = new double [length];

    // Initialize
    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k <= nz; k++) {
                dxF(i, j, k) = 0;
                dyF(i, j, k) = 0;
                dzF(i, j, k) = 0;
            }

    // dxF
    for (int i = delta; i <= nx + delta; i++)
        for (int j = 0; j <= ny + delta * 2; j++)
            for (int k = 0; k <= nz + delta * 2; k++) {
                S1(i, j, k) = 0;
                for (int di = -delta; di <= delta; di++)
                    S1(i, j, k) += _F(i + di, j, k) * (di * dx) * GaussianGradientCoefficient(sigma) * GaussianGradientExponent(di * dx, sigma);
            }
    for (int i = delta; i <= nx + delta; i++)
        for (int j = delta; j <= ny + delta; j++)
            for (int k = 0; k <= nz + delta * 2; k++) {
                S2(i, j, k) = 0;
                for (int dj = -delta; dj <= delta; dj++)
                    S2(i, j, k) += S1(i, j + dj, k) * GaussianGradientExponent(dj * dy, sigma);
            }
    for (int i = delta; i <= nx + delta; i++)
        for (int j = delta; j <= ny + delta; j++)
            for (int k = delta; k <= nz + delta; k++)
                for (int dk = -delta; dk <= delta; dk++)
                    dxF(i - delta, j - delta, k - delta) += S2(i, j, k + dk) * GaussianGradientExponent(dk * dz, sigma);
    // dyF
    for (int i = 0; i <= nx + delta * 2; i++)
        for (int j = delta; j <= ny + delta; j++)
            for (int k = 0; k <= nz + delta * 2; k++) {
                S1(i, j, k) = 0;
                for (int dj = -delta; dj <= delta; dj++)
                    S1(i, j, k) += _F(i, j + dj, k) * (dj * dy) * GaussianGradientCoefficient(sigma) * GaussianGradientExponent(dj * dy, sigma);
            }
    for (int i = delta; i <= nx + delta; i++)
        for (int j = delta; j <= ny + delta; j++)
            for (int k = 0; k <= nz + delta * 2; k++) {
                S2(i, j, k) = 0;
                for (int di = -delta; di <= delta; di++)
                    S2(i, j, k) += S1(i + di, j, k) * GaussianGradientExponent(di * dx, sigma);
            }
    for (int i = delta; i <= nx + delta; i++)
        for (int j = delta; j <= ny + delta; j++)
            for (int k = delta; k <= nz + delta; k++)
                for (int dk = -delta; dk <= delta; dk++)
                    dyF(i - delta, j - delta, k - delta) += S2(i, j, k + dk) * GaussianGradientExponent(dk * dz, sigma);
    // dxF
    for (int i = 0; i <= nx + delta * 2; i++)
        for (int j = 0; j <= ny + delta * 2; j++)
            for (int k = delta; k <= nz + delta; k++) {
                S1(i, j, k) = 0;
                for (int dk = -delta; dk <= delta; dk++)
                    S1(i, j, k) += _F(i, j, k + dk) * (dk * dz) * GaussianGradientCoefficient(sigma) * GaussianGradientExponent(dk * dz, sigma);
            }
    for (int i = delta; i <= nx + delta; i++)
        for (int j = 0; j <= ny + delta * 2; j++)
            for (int k = delta; k <= nz + delta; k++) {
                S2(i, j, k) = 0;
                for (int di = -delta; di <= delta; di++)
                    S2(i, j, k) += S1(i + di, j, k) * GaussianGradientExponent(di * dx, sigma);
            }
    for (int i = delta; i <= nx + delta; i++)
        for (int j = delta; j <= ny + delta; j++)
            for (int k = delta; k <= nz + delta; k++)
                for (int dj = -delta; dj <= delta; dj++)
                    dzF(i - delta, j - delta, k - delta) += S2(i, j + dj, k) * GaussianGradientExponent(dj * dy, sigma);

    delete [] stage1;
    delete [] stage2;
}

#define V(i, j, k) values[lcs::Code((i), (j), (k), ny + stencilSize - 1, nz + stencilSize - 1)]

int Crop(int lower, int value, int upper) {
    if (value < lower) return lower;
    if (value > upper) return upper;
    return value;
}

void DataExtension(double *values, int nx, int ny, int nz, int stencilSize) {
    int delta = stencilSize / 2;
    // Upper
    for (int i = 0; i <= nx + delta * 2; i++)
        for (int j = 0; j <= ny + delta * 2; j++)
            for (int k = 0; k <= nz + delta * 2; k++) {
                int _i = Crop(delta, i, nx + delta);
                int _j = Crop(delta, j, ny + delta);
                int _k = Crop(delta, k, nz + delta);
                V(i, j, k) = V(_i, _j, _k);
            } 
}

void GetSmoothedFTLE() {
    printf("Get the smoothed FTLE ...\n");

    int nx = configure.GetNx();
    int ny = configure.GetNy();
    int nz = configure.GetNz();
    int stencilSize = configure.GetStencilSize();
    int delta = stencilSize / 2;
    int length = (nx + 1) * (ny + 1) * (nz + 1);

    double dx = configure.GetDx();
    double dy = configure.GetDy();
    double dz = configure.GetDz();
    double sigma = configure.GetSigma();

    ftleValues = new double [length];
    N1 = new lcs::Vector [length];

    double *values = new double [(nx + stencilSize) * (ny + stencilSize) + (nz + stencilSize)];
    double *_dx_x = new double [length];
    double *_dy_x = new double [length];
    double *_dz_x = new double [length];
    double *_dx_y = new double [length];
    double *_dy_y = new double [length];
    double *_dz_y = new double [length];
    double *_dx_z = new double [length];
    double *_dy_z = new double [length];
    double *_dz_z = new double [length];

#define INIT_VALUES(dim) \
    for (int i = 0; i <= nx; i++) \
        for (int j = 0; j <= ny; j++) \
            for (int k = 0; k <= nz; k++) \
                values[lcs::Code(i + delta, j + delta, k + delta, ny + delta * 2, nz + delta * 2)] = P(i, j, k).Get##dim(); \
    DataExtension(values, nx, ny, nz, stencilSize);

    // dx_x, dy_x, dz_x
    INIT_VALUES(X)
    GetSmoothedDerivative(values, _dx_x, _dy_x, _dz_x, nx, ny, nz, dx, dy, dz, sigma, stencilSize);

    printf("Finished dx_x, dy_x, dz_x.\n");

    // dx_y, dy_y, dz_y
    INIT_VALUES(Y)
    GetSmoothedDerivative(values, _dx_y, _dy_y, _dz_y, nx, ny, nz, dx, dy, dz, sigma, stencilSize);

    printf("Finished dx_y, dy_y, dz_y.\n");

    // dx_z, dy_z, dz_z
    INIT_VALUES(Z)
    GetSmoothedDerivative(values, _dx_z, _dy_z, _dz_z, nx, ny, nz, dx, dy, dz, sigma, stencilSize);

    printf("Finished dx_z, dy_z, dz_z.\n");

    double minFTLE = 1e100, maxFTLE = -1e100;

    // Get FTLE and N1
    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k <= nz; k++) {
                double dx_x = _dx_x[lcs::Code(i, j, k, ny, nz)];
                double dy_x = _dy_x[lcs::Code(i, j, k, ny, nz)];
                double dz_x = _dz_x[lcs::Code(i, j, k, ny, nz)];
                double dx_y = _dx_y[lcs::Code(i, j, k, ny, nz)];
                double dy_y = _dy_y[lcs::Code(i, j, k, ny, nz)];
                double dz_y = _dz_y[lcs::Code(i, j, k, ny, nz)];
                double dx_z = _dx_z[lcs::Code(i, j, k, ny, nz)];
                double dy_z = _dy_z[lcs::Code(i, j, k, ny, nz)];
                double dz_z = _dz_z[lcs::Code(i, j, k, ny, nz)];
             
                double arr[] = {dx_x, dy_x, dz_x, dx_y, dy_y, dz_y, dx_z, dy_z, dz_z};

                lcs::Matrix F(arr, 3, 3);
                lcs::Matrix C = lcs::MatrixMatrixMultiplication(lcs::MatrixTranspose(F), F);

                double eigenvalues[3];
                double **eigenvectors = lcs::Matrix::CreateVTKMatrix(3, 3);
                double **CForVTK = C.MatrixForVTK();
                vtkMath::Jacobi(CForVTK, eigenvalues, eigenvectors);
                delete [] CForVTK;

                FTLE(i, j, k) = log(sqrt(eigenvalues[0]));
                N1(i, j, k) = lcs::Vector(eigenvectors[0][0], eigenvectors[1][0], eigenvectors[2][0]);

                lcs::Matrix::DisposeVTKMatrix(eigenvectors);

                minFTLE = std::min(minFTLE, FTLE(i, j, k));
                maxFTLE = std::max(maxFTLE, FTLE(i, j, k));
            }

    printf("minFTLE = %lf, maxFTLE = %lf\n", minFTLE, maxFTLE);
/*
    vtkSmartPointer<vtkImageData> image3D = vtkSmartPointer<vtkImageData>::New();
    image3D->SetExtent(0, nx, 0, ny, 0, nz);
    image3D->AllocateScalars(VTK_DOUBLE, 1);

    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k <= nz; k++)
                image3D->SetScalarComponentFromDouble(i, j, k, 0, FTLE(i, j, k));

    OutputImageData(image3D, nx, ny, nz, "SmoothedFTLEvalues.vtu");
*/

    vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
    data->SetNumberOfComponents(1);
    for (int k = 0; k <= nz; k++)
        for (int j = 0; j <= ny; j++)
            for (int i = 0; i <= nx; i++) {
                points->InsertNextPoint(i * dx, j * dy, k * dz);
                data->InsertNextTuple1(FTLE(i, j, k));
            }
    grid->SetDimensions(nx + 1, ny + 1, nz + 1);
    grid->SetPoints(points);
    grid->GetPointData()->SetScalars(data);
    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetInputData(grid);
    writer->SetFileName("SmoothedFTLEValue.vtu");
    writer->Write();

    delete [] values;
    delete [] _dx_x;
    delete [] _dy_x;
    delete [] _dz_x;
    delete [] _dx_y;
    delete [] _dy_y;
    delete [] _dz_y;
    delete [] _dx_z;
    delete [] _dy_z;
    delete [] _dz_z;

    printf("Done.\n\n");
}

void GetFakedFTLE() {
    printf("Get the faked FTLE ...\n");

    int nx = configure.GetNx();
    int ny = configure.GetNy();
    int nz = configure.GetNz();
    int length = (nx + 1) * (ny + 1) * (nz + 1);

    double dx = configure.GetDx();
    double dy = configure.GetDy();
    double dz = configure.GetDz();

    ftleValues = new double [length];
    N1 = new lcs::Vector [length];

    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k <= nz; k++) {
                double rule = nz / 2 + 0.5;
                if (k < rule) {
                    FTLE(i, j, k) = 300 + k - rule;
                    N1(i, j, k) = lcs::Vector(0, 0, 1);
                } else {
                    FTLE(i, j, k) = 300 + rule - k;
                    N1(i, j, k) = lcs::Vector(0, 0, -1);
                }
            }

    printf("Done.\n\n");
}

void InitializeValues(double *values, double *original, int nx, int ny, int nz, int stencilSize) {
    int delta = stencilSize / 2;

    printf("delta = %d, stencilSize = %d\n", delta, stencilSize);

    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k <= nz; k++)
                values[lcs::Code(i + delta, j + delta, k + delta, ny + delta * 2, nz + delta * 2)] = original[lcs::Code(i, j, k, ny, nz)];

    DataExtension(values, nx, ny, nz, stencilSize);
}

void GetSmoothedFTLEGradient() {
    printf("Get the smoothed FTLE gradient ...\n");

    int nx = configure.GetNx();
    int ny = configure.GetNy();
    int nz = configure.GetNz();
    int stencilSize = configure.GetStencilSize();
    int length = (nx + 1) * (ny + 1) * (nz + 1);

    double dx = configure.GetDx();
    double dy = configure.GetDy();
    double dz = configure.GetDz();
    double sigma = configure.GetSigma();

    ftleGradients = new lcs::Vector [length];

    double *values = new double [(nx + stencilSize) * (ny + stencilSize) * (nz + stencilSize)];
    double *_dx_ftle = new double [length];
    double *_dy_ftle = new double [length];
    double *_dz_ftle = new double [length];

    InitializeValues(values, ftleValues, nx, ny, nz, stencilSize);
    GetSmoothedDerivative(values, _dx_ftle, _dy_ftle, _dz_ftle, nx, ny, nz, dx, dy, dz, sigma, stencilSize);

    for (int i = 0; i < length; i++)
        ftleGradients[i] = lcs::Vector(_dx_ftle[i], _dy_ftle[i], _dz_ftle[i]);

    vtkSmartPointer<vtkImageData> image3D = vtkSmartPointer<vtkImageData>::New();
    image3D->SetExtent(0, nx, 0, ny, 0, nz);
    image3D->AllocateScalars(VTK_DOUBLE, 3);

    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k <= nz; k++) {
                image3D->SetScalarComponentFromDouble(i, j, k, 0, ftleGradients[lcs::Code(i, j, k, ny, nz)].GetX());
                image3D->SetScalarComponentFromDouble(i, j, k, 1, ftleGradients[lcs::Code(i, j, k, ny, nz)].GetY());
                image3D->SetScalarComponentFromDouble(i, j, k, 2, ftleGradients[lcs::Code(i, j, k, ny, nz)].GetZ());
            }

    OutputImageData(image3D, nx, ny, nz, "SmoothedFTLEGradients.vtu");

    delete [] values;
    delete [] _dx_ftle;
    delete [] _dy_ftle;
    delete [] _dz_ftle;

    printf("Done.\n\n");
}

#define __F(i, j, k) F[lcs::Code((i), (j), (k), ny, nz)]

void GetRawDerivative(double *F, double *dxF, double *dyF, double *dzF, int nx, int ny, int nz, double dx, double dy, double dz) {
    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k <= nz; k++) {
                if (!i)
                    dxF(i, j, k) = (__F(i + 1, j, k) - __F(i, j, k)) / dx;
                else if (i == nx)
                    dxF(i, j, k) = (__F(i, j, k) - __F(i - 1, j, k)) / dx;
                else
                    dxF(i, j, k) = (__F(i + 1, j, k) - __F(i - 1, j, k)) / (2 * dx);
                if (!j)
                    dyF(i, j, k) = (__F(i, j + 1, k) - __F(i, j, k)) / dy;
                else if (j == ny)
                    dyF(i, j, k) = (__F(i, j, k) - __F(i, j - 1, k)) / dy;
                else
                    dyF(i, j, k) = (__F(i, j + 1, k) - __F(i, j - 1, k)) / (2 * dy);
                if (!k)
                    dzF(i, j, k) = (__F(i, j, k + 1) - __F(i, j, k)) / dz;
                else if (k == nz)
                    dzF(i, j, k) = (__F(i, j, k) - __F(i, j, k - 1)) / dz;
                else
                    dzF(i, j, k) = (__F(i, j, k + 1) - __F(i, j, k - 1)) / ( 2 * dz);
            }
}

void GetRawFTLEGradient() {
    printf("Get the raw FTLE gradient ...\n");

    int nx = configure.GetNx();
    int ny = configure.GetNy();
    int nz = configure.GetNz();
    int length = (nx + 1) * (ny + 1) * (nz + 1);

    double dx = configure.GetDx();
    double dy = configure.GetDy();
    double dz = configure.GetDz();

    ftleGradients = new lcs::Vector [length];

    double *_dx_ftle = new double [length];
    double *_dy_ftle = new double [length];
    double *_dz_ftle = new double [length];

    GetRawDerivative(ftleValues, _dx_ftle, _dy_ftle, _dz_ftle, nx, ny, nz, dx, dy, dz);

    for (int i = 0; i < length; i++)
        ftleGradients[i] = lcs::Vector(_dx_ftle[i], _dy_ftle[i], _dz_ftle[i]);

    vtkSmartPointer<vtkImageData> image3D = vtkSmartPointer<vtkImageData>::New();
    image3D->SetExtent(0, nx, 0, ny, 0, nz);
    image3D->AllocateScalars(VTK_DOUBLE, 3);

    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k <= nz; k++) {
                image3D->SetScalarComponentFromDouble(i, j, k, 0, ftleGradients[lcs::Code(i, j, k, ny, nz)].GetX());
                image3D->SetScalarComponentFromDouble(i, j, k, 1, ftleGradients[lcs::Code(i, j, k, ny, nz)].GetY());
                image3D->SetScalarComponentFromDouble(i, j, k, 2, ftleGradients[lcs::Code(i, j, k, ny, nz)].GetZ());
            }

    OutputImageData(image3D, nx, ny, nz, "RawFTLEGradients.vtu");

    delete [] _dx_ftle;
    delete [] _dy_ftle;
    delete [] _dz_ftle;

    printf("Done.\n\n");
}

void PCATest() {
    printf("PCA Test ...\n");
    //double arr[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    //double arr[] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    //double arr[] = {3, 3, 3, -1, -1, -1, 2, 2, 2};
    //double arr[] = {1, 2, 3, 4, 8, 12, 0.5, 1, 1.5};
    //double arr[] = {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0};
    //double arr[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
    double arr[] = {1000, 999, 0, -999, -1000, 0};
    lcs::Matrix A(arr, 2, 3);
    lcs::Matrix B = lcs::PCA(A);
    B.Output();
    printf("Done.\n\n");
}

struct Triangle {
    lcs::Vector points[3];
    double ftleValues[3];
};

const int vertexList[8][3] = {
    {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}
};

const int edgeList[12][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0},
    {4, 5}, {5, 6}, {6, 7}, {7, 4},
    {0, 4}, {1, 5}, {2, 6}, {3, 7}
};

lcs::Vector EdgeInterpolation(const lcs::Vector &v1, const lcs::Vector &v2, double f1, double f2) {
    return (v2 * f1 - v1 * f2) / (f1 - f2);
}

double EdgeInterpolationScalar(double v1, double v2, double f1, double f2) {
    return (v2 * f1 - v1 * f2) / (f1 - f2);
}

void ClassifyCube(int a, int b, int c, int ny, int nz, double dx, double dy, double dz, double value[2][2][2], vtkSmartPointer<vtkPoints> surfacePoints, vtkSmartPointer<vtkDoubleArray> pointValues) {
    lcs::Vector corner(a * dx, b * dy, c * dz);

    lcs::Vector vertices[8];
    for (int i = 0; i < 8; i++)
        vertices[i] = corner + lcs::Vector(dx * vertexList[i][0],
                                           dy * vertexList[i][1],
                                           dz * vertexList[i][2]);

    double field[8];
    for (int i = 0; i < 8; i++)
        field[i] = value[vertexList[i][0]][vertexList[i][1]][vertexList[i][2]];

    int index = 0;
    for (int i = 0; i < 8; i++)
        index |= (field[i] < 0) << i;
        //index |= (FTLE(a + vertexList[i][0], b + vertexList[i][1], c + vertexList[i][2]) < 21) << i;

    lcs::Vector edgeVertices[12];
    lcs::Vector edgeVectors[12];
    double ftles[12];
    for (int i = 0; i < 12; i++) {
        int vtx1 = edgeList[i][0];
        int vtx2 = edgeList[i][1];
        edgeVertices[i] = EdgeInterpolation(vertices[vtx1], vertices[vtx2], field[vtx1], field[vtx2]);
        //double ftle1 = FTLE(a + vertexList[vtx1][0], b + vertexList[vtx1][1], c + vertexList[vtx1][2]) - 21;
        //double ftle2 = FTLE(a + vertexList[vtx2][0], b + vertexList[vtx2][1], c + vertexList[vtx2][2]) - 21;
        //ftles[i] = EdgeInterpolationScalar(ftle1, ftle2, field[vtx1], field[vtx2]);
        //edgeVertices[i] = EdgeInterpolation(vertices[vtx1], vertices[vtx2], ftle1, ftle2);

        //lcs::Vector norm1 = N1(a + vertexList[vtx1][0], b + vertexList[vtx1][1], c + vertexList[vtx1][2]);
        //lcs::Vector norm2 = N1(a + vertexList[vtx2][0], b + vertexList[vtx2][1], c + vertexList[vtx2][2]);
        lcs::Vector grad1 = ftleGradients[lcs::Code(a + vertexList[vtx1][0], b + vertexList[vtx1][1], c + vertexList[vtx1][2], ny, nz)];
        lcs::Vector grad2 = ftleGradients[lcs::Code(a + vertexList[vtx2][0], b + vertexList[vtx2][1], c + vertexList[vtx2][2], ny, nz)];
        //edgeNormals[i] = EdgeInterpolation(grad1, grad2, ftle1, ftle2);
        //ftles[i] = 21;
        edgeVectors[i] = EdgeInterpolation(grad1, grad2, field[vtx1], field[vtx2]);
    }

    int numOfVertices = numVertsTable[index];

    for (int i = 0; i < numOfVertices; i += 3) {
        //Triangle triangle;
        for (int j = 0; j < 3; j++) {
            int edgeIdx = triTable[index][i + j];

            surfacePoints->InsertNextPoint(edgeVertices[edgeIdx].GetX(),
                                           edgeVertices[edgeIdx].GetY(),
                                           edgeVertices[edgeIdx].GetZ());
            pointValues->InsertNextTuple3(edgeVectors[edgeIdx].GetX(), edgeVectors[edgeIdx].GetY(), edgeVectors[edgeIdx].GetZ());
            //pointValues->InsertNextTuple1(ftles[edgeIdx]);
            //triangle.points[j] = edgeVertices[edgeIdx];
            //triangle.ftleValues[j] = ftleValues[edgeIdx];
        }
        //triangles.push_back(triangle);
    }
}

void SurfaceExtraction() {
    printf("Extract the surface ...\n");

    int nx = configure.GetNx();
    int ny = configure.GetNy();
    int nz = configure.GetNz();
    int length = (nx + 1) * (ny + 1) * (nz + 1);

    double dx = configure.GetDx();
    double dy = configure.GetDy();
    double dz = configure.GetDz();

    double arr[24];
    double scalar[2][2][2];

    vtkSmartPointer<vtkPoints> surfacePoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> pointValues = vtkSmartPointer<vtkDoubleArray>::New();
    //pointValues->SetNumberOfComponents(1);
    pointValues->SetNumberOfComponents(3);

    int tot = 0, blk = 0;

    printf("Check point 1\n");

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                // if (++tot % 1000 == 0) printf("tot = %d, blk = %d, blk / tot = %lf, numOfTri = %lf\n", tot, blk, (double)blk / tot, (double)surfacePoints->GetNumberOfPoints() / 3);

                /// DEBUG ///
                if (i == 10 && j == 20 && k == 30) {
                  for (int di = 0; di < 2; di++) {
                    for (int dj = 0; dj < 2; dj++) {
                      for (int dk = 0; dk < 2; dk++) {
                        printf("v: %lf\n", FTLE(i + di, j + dj, k + dk));
                      }
                    }
                  }
                }

                int cnt = 0;
                for (int di = 0; di <= 1; di++)
                    for (int dj = 0; dj <= 1; dj++)
                        for (int dk = 0; dk <= 1; dk++) {
                            arr[cnt++] = N1(i + di, j + dj, k + dk).GetX();
                            arr[cnt++] = N1(i + di, j + dj, k + dk).GetY();
                            arr[cnt++] = N1(i + di, j + dj, k + dk).GetZ();
                        }
/*
                arr[cnt++] = 0;
                arr[cnt++] = 0;
                arr[cnt++] = 0;
                lcs::Matrix A(arr, 9, 3);
*/
                lcs::Matrix A(arr, 8, 3);
                A = lcs::PCA(A);
                //lcs::Vector rule(A.Element(0, 0), A.Element(0, 1), A.Element(0, 2));
                lcs::Vector rule(A.Element(0, 0), A.Element(1, 0), A.Element(2, 0));

                /// DEBUG ///

                int lower = 0;

                for (int di = 0; di <= 1; di++)
                    for (int dj = 0; dj <= 1; dj++)
                        for (int dk = 0; dk <= 1; dk++) {
                            lcs::Vector e = N1(i + di, j + dj, k + dk);
                            if (lcs::Dot(e, rule) < 0) e = e * -1;
                            scalar[di][dj][dk] = lcs::Dot(e, ftleGradients[lcs::Code(i + di, j + dj, k + dk, ny, nz)]);

                            /// DEBUG ///
                            if (i == 10 && j == 20 && k == 30) {
                                printf("e: %lf %lf %lf\n", e.GetX(), e.GetY(), e.GetZ());
                                lcs::Vector g = ftleGradients[lcs::Code(i + di, j + dj, k + dk, ny, nz)];
                                printf("g: %lf %lf %lf\n", g.GetX(), g.GetY(), g.GetZ());
                                printf("scalar: %lf\n", scalar[di][dj][dk]);
                                printf("\n");
                            }

                            lower += FTLE(i + di, j + dj, k + dk) < 21.0;
                        }

                /// primary filtering
                // if (lower < 8) {
                    ClassifyCube(i, j, k, ny, nz, dx, dy, dz, scalar, surfacePoints, pointValues);
                    blk++;
                // }
            }

    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    for (vtkIdType i = 0; i < surfacePoints->GetNumberOfPoints() / 3; i++) {
        cells->InsertNextCell(3);
        cells->InsertCellPoint(i * 3);
        cells->InsertCellPoint(i * 3 + 1);
        cells->InsertCellPoint(i * 3 + 2);
    }
    vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
    poly->SetPoints(surfacePoints);
    poly->SetPolys(cells);
    //poly->GetPointData()->SetScalars(pointValues);
    poly->GetPointData()->SetVectors(pointValues);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetInputData(poly);
    writer->SetFileName("surface.vtu");
    writer->Write();

    printf("Done.\n\n");
}

int main() {
    PCATest();
    ReadConfiguration();
    ReadLastPositions();
    GetRawFTLE(); // It is just for comparison to smoothed FTLE.
    //GetSmoothedFTLE(); // totally separate from GetRawFTLE()
    //GetFakedFTLE();
    //GetSmoothedFTLEGradient();
    GetRawFTLEGradient();
    SurfaceExtraction();
    return 0;
}
