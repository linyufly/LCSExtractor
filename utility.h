/**********************************************
File		:	utility.h
Author		:	Mingcheng Chen
Last Update	:	January 1st, 2014
***********************************************/

#ifndef __UTILITY_H
#define __UTILITY_H

namespace lcs {

void Error(const char *format, ...);
double ParseDouble(const char *str);
int ParseInt(const char *str);

class Configuration {
private:
    double dx;
    double dy;
    double dz;
    double sigma; // for smoothing
    int nx;
    int ny;
    int nz;

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

    int GetNx() const {
        return nx;
    }

    int GetNy() const {
        return ny;
    }

    int GetNz() const {
        return nz;
    }
};

}

#endif
