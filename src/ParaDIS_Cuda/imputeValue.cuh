#ifndef IMPUTEVALUE_H
#define IMPUTEVALUE_H

#include <malloc.h>
#include <omp.h>
#include <limits.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "DTW.h"
#include "Z_norm.h"
#include "debugger.h"
#include "profiler.h"

itemType impute(itemType* S, itemType** R, int h, itemType *SetC_gpu, itemType *Q_gpu, itemType *LB_gpu);
void verticalOverlap(int h, bool corr);
void FillData(itemType* S, itemType** R, int h);
void CslcZnorm(int h);
int CalcLB(int h, itemType *SetC_gpu, itemType *Q_gpu, itemType *LB_gpu);
__global__ void LB_calc_kernel(itemType *Q, itemType *SetC, itemType *LB, int N, int idx_max);
__host__ __device__ itemType LbKim(itemType* Q, itemType* C);
__host__ __device__ itemType LbKeogh(itemType* Q, itemType* C, int r);

#endif // IMPUTEVALUE_H



