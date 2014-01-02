/**********************************************
File		:	utility.cpp
Author		:	Mingcheng Chen
Last Update	:	January 1st, 2014
***********************************************/

#include <cstdio>
#include <cstdlib>
#include <cstdarg>

#include "utility.h"

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

void lcs::Configuration::LoadFile(const char *fileName) {
    FILE *fin;
    fin = fopen(fileName, "r");
    
    fclose(fin);
}
