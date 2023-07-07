#ifndef READDATA_H
#define READDATA_H

#include <assert.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"



void read_timeseries(const char *filename, itemType *S, itemType**R);
const char* getelem(char* line, int num);

#endif // READDATA_H
