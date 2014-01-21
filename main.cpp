/**********************************************
File		:	main.cpp
Author		:	Mingcheng Chen
Last Update	:	January 14th, 2014
***********************************************/

#include "lcsGeometry.h"
#include "utility.h"
#include "algebra.h"

#include <cmath>

#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>

lcs::Configuration configure;
lcs::Vector *lastPositions;
double *ftleValues;
lcs::Vector *N1;

void ReadConfiguration() {
    configure.LoadFile("configuration.conf");
    printf("Finished reading the configuration file.\n");
    printf("\n");
}

void ReadLastPositions() {
    printf("Start reading the data file ... \n");

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
            }

    vtkSmartPointer<vtkImageData> image3D = vtkSmartPointer<vtkImageData>::New();
    image3D->SetExtent(0, nx, 0, ny, 0, nz);
    image3D->AllocateScalars(VTK_DOUBLE, 1);

    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k <= nz; k++)
                image3D->SetScalarComponentFromDouble(i, j, k, 0, FTLE(i, j, k));

    OutputImageData(image3D, nx, ny, nz, "SmoothedFTLEvalues.vtu");

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

int main() {
    ReadConfiguration();
    ReadLastPositions();
    GetRawFTLE();
    GetSmoothedFTLE(); // Totally separate from GetRawFTLE()
    return 0;
}
