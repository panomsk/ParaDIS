/*+++++++++++++++++++++++++++++++++

Project: ParaDIS (Parallel Algorithm for Imputation of Missing Values in Streaming Time Series)

Source file: Z_norm.h

Purpose: Z-normalization of subsequence

Author(s): Andrey Poluyanov (andrey.poluyanov@gmail.com) and Mikhail Zymbler (mzym@susu.ru)

+++++++++++++++++++++++++++++++++*/

#ifndef Z_NORM_H
#define Z_NORM_H

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "config.h"
#include "ReadData.h"

void Z_norm(itemType* query, int lenQuery);


#endif // Z_NORM_H

