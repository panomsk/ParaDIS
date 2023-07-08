/*+++++++++++++++++++++++++++++++++

Project: ParaDIS (Parallel Algorithm for Imputation of Missing Values in Streaming Time Series)

Source file: imputeValue.h

Purpose: Imputation one missing point in time series 

Author(s): Andrey Poluyanov (andrey.poluyanov@gmail.com) and Mikhail Zymbler (mzym@susu.ru)

+++++++++++++++++++++++++++++++++*/

#ifndef IMPUTEVALUE_H
#define IMPUTEVALUE_H

#include <malloc.h>
#include <omp.h>
#include <limits.h>

#include "config.h"
#include "DTW.h"
#include "lowerBounds.h"
#include "Z_norm.h"
#include "debugger.h"
#include "profiler.h"

itemType impute(itemType* S, itemType** R, int h);
void verticalOverlap(int h, bool corr);
void FillData(itemType* S, itemType** R, int h);
void CslcZnorm(int h);
int CalcLB(int h);

#endif // IMPUTEVALUE_H



