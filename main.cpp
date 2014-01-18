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
    vtkSmartPointer<vtkXMLImageDataWriter> writer =
    vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(fileName);

#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(data->GetProducerPort());
#else
    writer->SetInputData(data);
#endif

    writer->Write();
}

void GetRawFTLE() {
    int nx = configure.GetNx();
    int ny = configure.GetNy();
    int nz = configure.GetNz();
    int length = (nx + 1) * (ny + 1) * (nz + 1);

    double dx = configure.GetDx();
    double dy = configure.GetDy();
    double dz = configure.GetDz();

    ftleValues = new double [length];
    N1 = new lcs::Vector [length];

    int i, j, k;
    for (i = 0; i <= nx; i++)
        for (j = 0; j <= ny; j++)
            for (k = 0; k <= nz; k++) {
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

                C.Output();

                double eigenvalues[3];
                double **eigenvectors = lcs::Matrix::CreateVTKMatrix(3, 3);
                double **CForVTK = C.MatrixForVTK();
                vtkMath::Jacobi(CForVTK, eigenvalues, eigenvectors);
                delete [] CForVTK;

                FTLE(i, j, k) = log(sqrt(eigenvalues[0]));
                printf("%d, %d, %d: %lf\n", i, j, k, FTLE(i, j, k));
                N1(i, j, k) = lcs::Vector(eigenvectors[0][0], eigenvectors[1][0], eigenvectors[2][0]);

                lcs::Matrix::DisposeVTKMatrix(eigenvectors);
            }

    vtkSmartPointer<vtkImageData> image3D = vtkSmartPointer<vtkImageData>::New();
    image3D->SetExtent(0, nx, 0, ny, 0, nz);
    image3D->AllocateScalars(VTK_DOUBLE,1);

    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k <= nz; k++)
                image3D->SetScalarComponentFromDouble(i, j, k, 0, FTLE(i, j, k));

    OutputImageData(image3D, nx, ny, nz, "FTLEvalues.vtu");
}

int main() {
    ReadConfiguration();
    ReadLastPositions();
    GetRawFTLE();
    return 0;
}
