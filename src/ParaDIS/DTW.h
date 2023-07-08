/*+++++++++++++++++++++++++++++++++

Project: ParaDIS (Parallel Algorithm for Imputation of Missing Values in Streaming Time Series)

Source file: DTW.h

Purpose: DTW calculation for two subsequences

Author(s): Andrey Poluyanov (andrey.poluyanov@gmail.com) and Mikhail Zymbler (mzym@susu.ru)

+++++++++++++++++++++++++++++++++*/



#ifndef DTW_H
#define DTW_H

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "config.h"

itemType DTW(itemType* Q, itemType* C, int r);

#endif // DTW_H