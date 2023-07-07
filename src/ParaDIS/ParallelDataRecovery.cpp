#include <time.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <omp.h>
#include <math.h>
#include <iostream>

#include "ReadData.h"
#include "Z_norm.h"
#include "config.h"
#include "DTW.h"
#include "imputeValue.h"
#include "debugger.h"
#include "profiler.h"

using namespace std;

char* data_fname;
int num_of_threads;
int DTW_count=0;
int LB_count=0;

int main(int argc, char * argv[])
{
	START("ParallelDataRecovery");
	PRF_INIT;
	
	struct timespec begin, end;
	clock_gettime(CLOCK_MONOTONIC_RAW, &begin);
	double start = omp_get_wtime();  
	
	data_fname = argv[1];
	assert(data_fname != NULL); 

	assert(argv[2]!=NULL);
    num_of_threads = atoi(argv[2]);
    assert(num_of_threads >= 1 && num_of_threads <= omp_get_max_threads());   
	
	itemType* S = (itemType*)malloc((L_par + numimputevalue) * sizeof(itemType));

	itemType** R = (itemType**)calloc(d_par, sizeof(itemType*));
	for (int i = 0; i < d_par; i++) {
		R[i] = (itemType*)malloc((L_par + numimputevalue) * sizeof(itemType));
	}

    PRF_START(start1);
	read_timeseries(data_fname, S, R);
	PRF_FINISH(finish1);
	
	itemType S_predicted[numimputevalue];
	itemType S_init[numimputevalue];

	for (int i = 0; i < numimputevalue; i++)
	{
		itemType s = impute(S, R, i);
		S_init[i] = S[L_par + i];
		S_predicted[i] = s;

		if(cum_error==true)
		{
		    S[L_par + i] = s;
		}
	}

	printf("RESULTS for csv export\n");
	printf("S_init;S_predicted\n");
	for (int i = 0; i < numimputevalue; i++)
	{
		printf("%lf;%lf\n", S_init[i], S_predicted[i]);
	}	
		
	clock_gettime(CLOCK_MONOTONIC_RAW, &end);
	double stop = omp_get_wtime();

	#ifdef PROFILE
		printf("Params_IDK:\n");
		printf("data_fname;cum_error;L_par;numimputevalue;init_bsf;cut_off_neighbors;num_lb;d_par;l_par;k_par;r_par;use_corr\n");
		printf("%s;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d\n",data_fname,cum_error,L_par,numimputevalue,init_bsf,cut_off_neighbors,num_lb,d_par,l_par,k_par,r_par,use_corr);

		printf("Time for reading data from file = %.2lf%%\n", (finish1 - start1) / (stop - start) * 100);
		printf("Time for FillData = %.2lf%%\n", (finish2 - start2) / (stop - start) * 100);
		printf("Time for CslcZnorm = %.2lf%%\n", (finish3 - start3) / (stop - start) * 100);
		printf("Time for LB = %.2lf%%\n", (finish4 - start4) / (stop - start) * 100);
		printf("Time for BitMap = %.2lf%%\n", (finish5 - start5) / (stop - start) * 100);
		printf("Time for DTW = %.2lf%%\n", (finish6 - start6) / (stop - start) * 100);
		printf("Time for rating and S-value = %.2lf%%\n", (finish7 - start7) / (stop - start) * 100);
		
		double check_time = 0;
		check_time += finish1 - start1;
		check_time += finish2 - start2;
		check_time += finish3 - start3;
		check_time += finish4 - start4;
		check_time += finish5 - start5;
		check_time += finish6 - start6;
		check_time += finish7 - start7;
		printf("Sum of previous steps = %.2lf%%\n", (check_time) / (stop - start) * 100);

		printf("Total Runtime counted by omp_get_wtime %.2lf sec\n", stop - start);
		printf("Total Runtime counted by clock_gettime %.2lf sec\n", (end.tv_nsec - begin.tv_nsec) / 1000000000.0 + (end.tv_sec  - begin.tv_sec));

		check_time -= finish1 - start1;
		check_time -= finish2 - start2;
		printf("Runtime without Read and FillData %.2lf sec\n", check_time);

		double MAE, MSE, RMSE;
		MAE=MSE=RMSE=0;
		for (int i = 0; i < numimputevalue; i++)
		{
			MAE += fabs(S_init[i]-S_predicted[i]);
			MSE += pow(S_init[i]-S_predicted[i], 2);
		}
		MAE = MAE/numimputevalue;
		MSE = MSE/numimputevalue;
		RMSE = sqrt(MSE);
		printf("MAE = %.5lf\n", MAE);
		printf("MSE = %.5lf\n", MSE);
		printf("RMSE = %.5lf\n", RMSE);
	#endif

	#ifdef DEBUG
		if(LB_count>0)
		{
			printf("LB count of calculating = %d\n", LB_count);
			printf("DTW count of calculating = %d\n", DTW_count);
			printf("DTW part of counting %.2lf %%\n", double(DTW_count)/LB_count*100);
		}
	#endif

	FINISH("ParallelDataRecovery");
	printf("\n\n");
	return 0;
}

