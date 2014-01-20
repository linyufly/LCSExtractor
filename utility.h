/**********************************************
File		:	utility.h
Author		:	Mingcheng Chen
Last Update	:	January 14th, 2014
***********************************************/

#ifndef __UTILITY_H
#define __UTILITY_H

#include <cstdio>

#include <string>

namespace lcs {

void Error(const char *format, ...);
double ParseDouble(const char *str);
int ParseInt(const char *str);
void ConsumeChar(char goal, FILE *fin);
double GetCurrentTimeInSeconds();

int Code(int x, int y, int z, int ny, int nz);
void Decode(int &x, int &y, int &z, int code, int ny, int nz);

class Configuration {
private:
    double dx;
    double dy;
    double dz;
    double sigma; // for smoothing
    int nx;
    int ny;
    int nz;
    int stencilSize;
    std::string dataFile;

public:
    void LoadFile(const char *fileName);

    double GetDx() const {
        return dx;
    }

    double GetDy() const {
        return dy;
    }

    double GetDz() const {
        return dz;
    }

    double GetSigma() const {
        return sigma;
    }

    std::string GetDataFile() const {
        return dataFile;
    }

    int GetNx() const {
        return nx;
    }

    int GetNy() const {
        return ny;
    }

    int GetNz() const {
        return nz;
    }

    int GetStencilSize() const {
        return stencilSize;
    }
};

}

#endif
