/*+++++++++++++++++++++++++++++++++

Project: ParaDIS (Parallel Algorithm for Imputation of Missing Values in Streaming Time Series)

Source file: DTW.cpp

Purpose: DTW calculation for two subsequences

Author(s): Andrey Poluyanov (andrey.poluyanov@gmail.com) and Mikhail Zymbler (mzym@susu.ru)

+++++++++++++++++++++++++++++++++*/

#include "DTW.h"

itemType DTW(itemType* Q, itemType* C, int r)
{
    itemType *cost;
    itemType *cost_prev;
    itemType *cost_tmp;

    cost = (itemType*)malloc(sizeof(itemType)*(l_par));
    cost_prev = (itemType*)malloc(sizeof(itemType)*(l_par));

    int i, j;

    for (i = 0; i <= (l_par-1); i++)
    {
        //cost[i] = inf_val;
        cost_prev[i] = inf_val;
    }

    cost_prev[0] = dist(Q[0],C[0]);
    
    for (j = 1; j <= min(l_par-1,r); j++)
    {
        cost_prev[j] = cost_prev[j-1] + dist(Q[0],C[j]);
    }

 
    for (i = 1; i<=(l_par-1); i++)
    {
        for (j = 0; j <= (l_par-1); j++)
        {
            cost[j] = inf_val;
        }

        for (j = max(0,i-r); j <= min(l_par-1,i+r); j++)
        {
            if ((j-1) < 0)
                cost[j] = cost_prev[j];
            else
                cost[j] = min(cost_prev[j],cost_prev[j-1]);
        }

        for (j = max(0,i-r); j <= min(l_par-1,i+r); j++)
        {
            if ((j-1) < 0)
                cost[j] = dist(Q[i],C[j]) + cost[j];
            else
                cost[j] = dist(Q[i],C[j]) + min(cost[j],cost[j-1]);            
        }
        
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;

    }

    itemType final_dtw = cost_prev[l_par-1];

    free(cost);
    free(cost_prev);

    return final_dtw;
} 