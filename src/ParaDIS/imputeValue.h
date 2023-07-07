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



