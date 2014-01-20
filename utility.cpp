/**********************************************
File		:	utility.cpp
Author		:	Mingcheng Chen
Last Update	:	January 14th, 2014
***********************************************/

#include "utility.h"

#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cctype>
#include <cstring>
#include <sys/time.h>

int lcs::Code(int x, int y, int z, int ny, int nz) {
    return (x * (ny + 1) + y) * (nz + 1) + z;
}

void lcs::Decode(int &x, int &y, int &z, int code, int ny, int nz) {
    z = code % (nz + 1);
    code /= (nz + 1);
    y = code % (ny + 1);
    x = code / (ny + 1);
}

void lcs::Error(const char *format, ...) {
    char buffer[256];
    va_list args;
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);
    printf("Error: %s\n", buffer);
    exit(-1);
}

double lcs::ParseDouble(const char *str) {
    double res;
    sscanf(str, "%lf", &res);
    return res;
}

int lcs::ParseInt(const char *str) {
    int res;
    sscanf(str, "%d", &res);
    return res;
}

void lcs::ConsumeChar(char goal, FILE *fin) {
    char ch = fgetc(fin);
    for (; ch != goal && ch != EOF; ch = fgetc(fin));
    if (ch == EOF) lcs::Error("The configuration file is defective.");
}

double lcs::GetCurrentTimeInSeconds() {
    timeval currTime;
    gettimeofday(&currTime, 0);
    return currTime.tv_sec + currTime.tv_usec * 1e-6;
}

#define QUOTE(str) #str
#define STR(str) QUOTE(str)
#define READ_INT(var) \
    if (!strcmp(name, STR(var))) { \
        printf(STR(read var ...)); \
        int value; \
        if (fscanf(fin, "%d", &value) != 1) lcs::Error(STR(Fail to read var)); \
        this->var = value; \
        printf(STR(Done. var = %d\n), var); \
        continue; \
    } 
#define READ_DOUBLE(var) \
    if (!strcmp(name, STR(var))) { \
        printf(STR(read var ...)); \
        double value; \
        if (fscanf(fin, "%lf", &value) != 1) lcs::Error(STR(Fail to read var)); \
        this->var = value; \
        printf(STR(Done. var = %lf\n), var); \
        continue; \
    } 
#define READ_STRING(var) \
    if (!strcmp(name, STR(var))) { \
        printf(STR(read var ...)); \
        lcs::ConsumeChar('\"', fin); \
        this->var = ""; \
        while (1) { \
            ch = fgetc(fin); \
            if (ch == EOF) lcs::Error("The configuration file is defective."); \
            if (ch == '\"') break; \
            this->var += ch; \
        } \
        printf(STR(Done. var = %s\n), var.c_str()); \
    }

void lcs::Configuration::LoadFile(const char *fileName) {
    FILE *fin;
    fin = fopen(fileName, "r");
 
    if (fin == NULL)
        lcs::Error("Configuration file \"%s\" cannot be opened.", fileName);

    char ch;
    while ((ch = fgetc(fin)) != EOF) {
        // Skip the comments
        if (ch == '#') {
            for (; ch != '\n' && ch != '\r' && ch != EOF; ch = fgetc(fin));
            continue;
        }

        // Get the variable name
        if (isalpha(ch)) {
            char name[50];
            int len = 0;
            for (; isalpha(ch); name[len++] = ch, ch = fgetc(fin));
            name[len] = 0;

            // Skip the equation mark
            lcs::ConsumeChar('=', fin);

            READ_DOUBLE(dx)
            READ_DOUBLE(dy)
            READ_DOUBLE(dz)
            READ_DOUBLE(sigma)
            READ_INT(nx)
            READ_INT(ny)
            READ_INT(nz)
            READ_INT(stencilSize)
            READ_STRING(dataFile)

            if (stencilSize % 2 == 0)
                lcs::Error("Stencil size is even.");
        }
    }

    fclose(fin);
}
