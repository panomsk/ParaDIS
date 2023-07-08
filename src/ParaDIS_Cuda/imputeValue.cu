/*+++++++++++++++++++++++++++++++++

Project: ParaDIS (Parallel Algorithm for Imputation of Missing Values in Streaming Time Series)

Source file: imputeValue.cu

Purpose: Imputation one missing point in time series 

Author(s): Andrey Poluyanov (andrey.poluyanov@gmail.com) and Mikhail Zymbler (mzym@susu.ru)

+++++++++++++++++++++++++++++++++*/

#include "imputeValue.cuh"

itemType Q[d_par][l_par];
itemType SetC[d_par][subseq_count + numimputevalue][l_par];
itemType LB[d_par][subseq_count + numimputevalue][num_lb];
itemType bsf[d_par];
itemType N[d_par][subseq_count + numimputevalue];
int CandIndex[subseq_count + numimputevalue];
int ind[k_par];
itemType RANK[subseq_count + numimputevalue];
bool BitMap[d_par][subseq_count + numimputevalue];
extern int num_of_threads;
extern int DTW_count;
extern int LB_count;

struct maxRank{
	itemType Value;
	int Indax;
};
#pragma omp declare reduction(max : struct maxRank : \
	omp_out.Value = omp_in.Value < omp_out.Value ? omp_out.Value : omp_in.Value, \
	omp_out.Indax = omp_in.Value < omp_out.Value ? omp_out.Indax : omp_in.Indax ) \
	initializer( omp_priv = { 0, 0 } )

#pragma omp declare reduction(min : struct maxRank : \
	omp_out.Value = omp_in.Value > omp_out.Value ? omp_out.Value : omp_in.Value, \
	omp_out.Indax = omp_in.Value > omp_out.Value ? omp_out.Indax : omp_in.Indax ) \
	initializer( omp_priv = { inf_val, 0 } )



itemType impute(itemType* S, itemType** R, int h, itemType *SetC_gpu, itemType *Q_gpu, itemType *LB_gpu)
{

    itemType s = 0;
	PRF_START(start2);
	FillData(S, R, h);
	PRF_FINISH(finish2);
	PRF_START(start3);
	CslcZnorm(h);
	PRF_FINISH(finish3);
	CalcLB(h, SetC_gpu, Q_gpu, LB_gpu);
	
	PRF_START(start7);
	verticalOverlap(h, use_corr);
	int i, j = 0;
	int rank_counts = 0;
	while (j < k_par)
	{
		struct maxRank maxRANK = {0, 0};
		#pragma omp parallel for num_threads(num_of_threads) reduction (max: maxRANK)
		for (i = 0; i <= subseq_count + h; i++)
		{
			if (maxRANK.Value < RANK[i])
			{ 
				maxRANK.Value = RANK[i];
				maxRANK.Indax = i;
			}
		}

		if(maxRANK.Value>0)
		{
			ind[j] = maxRANK.Indax;
			int t;
			for (t = max(0, (ind[j] - l_par)); t < min((ind[j] + l_par),(subseq_count + h + 1)); t++)
			{
				RANK[t] = 0;
			}
			rank_counts+=1;
		}
		else break;
		j++;
	}

	assert(rank_counts>0);

	#ifdef DEBUG
		DTW_count += d_par;
		LB_count += d_par*(subseq_count + h + 1);	
		for (j = 0; j < d_par; j++)
		{
			for (int i = 0; i <= subseq_count + h; i++)
			{
				if (N[j][i] < inf_val)
				{
					DTW_count++;
				}
			}
		}
	#endif

	int n = 0;
	for (j = 0; j < k_par; j++) {
		if (ind[j] != -1) 
		{
			s = s + S[ind[j] + l_par -1];
			n++;
		}
	}
	assert(n>0);
	
	PRF_FINISH(finish7);
	return s/n;
}

void FillData(itemType* S, itemType** R, int h)
{
	int i, j, q;
	for (j = 0; j < d_par; j++)
	{
		for (i = 0; i < l_par; i++)
			Q[j][i] = R[j][L_par - l_par + h + 1 + i];
	}
	for (i = 0; i < d_par; i++)
	{
		for (j = 0; j <= subseq_count + h; j++)
		{
			for (q = 0; q < l_par; q++)
				SetC[i][j][q] = R[i][j + q];
		}
	}
	for (i = 0; i < d_par; i++)
	{
		for (j = 0; j <= subseq_count + h; j++)
			N[i][j] = inf_val;
	}

	for (i = 0; i < k_par; i++)
		ind[i] = -1;
	
	for (j = 0; j <= subseq_count + h; j++)
	{
		RANK[j] = 0;
	}

	for (i = 0; i < d_par; i++)
	{
		for (j = 0; j <= subseq_count + h; j++)
			BitMap[i][j] = true;
	}
}

void CslcZnorm(int h)
{
	int i ,j;
	#pragma omp parallel for num_threads(num_of_threads) collapse(2) 
	for (i = 0; i < d_par; i++)
	{
		for (j = 0; j <= subseq_count + h; j++) 
			Z_norm(SetC[i][j], l_par);
	}
	#pragma omp parallel for num_threads(num_of_threads)
	for (i = 0; i < d_par; i++)
	{
		Z_norm(Q[i], l_par);
	}
}

int CalcLB(int h, itemType *SetC_gpu, itemType *Q_gpu, itemType *LB_gpu)
{
	int i,j,s,cnt,left,right;
	itemType cur_dist;
	PRF_START(start4);
	
	
	if (use_GPU==true)
	{
		int N_subs = subseq_count + numimputevalue;
		
		cudaMemcpy(SetC_gpu, SetC, d_par*(subseq_count + numimputevalue)*l_par*sizeof(itemType), cudaMemcpyHostToDevice);
		cudaMemcpy(Q_gpu, Q, d_par*l_par*sizeof(itemType), cudaMemcpyHostToDevice);
		
		dim3 grid, block;
		int blockSize = 512;
		
		grid.x = (N_subs + blockSize - 1) / blockSize;  grid.y = d_par;
		block.x = blockSize; block.y = 1;
		
		LB_calc_kernel<<<grid, block>>>(Q_gpu, SetC_gpu, LB_gpu, N_subs, (subseq_count + h + 1)); 
		cudaMemcpy(LB, LB_gpu, d_par*(subseq_count + numimputevalue)*num_lb*sizeof(itemType), cudaMemcpyDeviceToHost);
	}
	else
	{		
	
		//CPU
		#pragma omp parallel for num_threads(num_of_threads) collapse(2)
		for (i = 0; i < d_par; i++)
		{
			for (j = 0; j <= subseq_count + h; j++)
			{
				LB[i][j][0] = LbKim(Q[i], SetC[i][j]);
				LB[i][j][1] = LbKeogh(Q[i], SetC[i][j], r_par);
				LB[i][j][2] = LbKeogh(SetC[i][j], Q[i], r_par);
				LB[i][j][3] = max(LB[i][j][1], LB[i][j][2]);
			}
		}
	
	}
	
	PRF_FINISH(finish4);


	for (int i = 0; i < d_par; i++)
	{
		bsf[i] = inf_val;
		if (init_bsf==0)
		{
			bsf[i] = DTW(Q[i], SetC[i][(L_par - (3 * l_par) + h)], r_par) + itemtype_epsilon;
		}
		else if (init_bsf==1)
		{
			bsf[i] = DTW(Q[i], SetC[i][subseq_count + h], r_par) + itemtype_epsilon;
		}
		else if (init_bsf==2)
		{
			struct maxRank minLB = {inf_val, 0};
			#pragma omp parallel for num_threads(num_of_threads) reduction (min: minLB)
			for (j = 0; j <= subseq_count + h; j++)
			{
				if (minLB.Value > LB[i][j][3])
				{ 
					minLB.Value = LB[i][j][3];
					minLB.Indax = j;
				}
			}
			bsf[i] = DTW(Q[i], SetC[i][minLB.Indax], r_par) + itemtype_epsilon;
		}
	}
		
	for (i = 0; i < d_par; i++)
	{
		while(true)
		{
			PRF_START(start5);
			#pragma omp parallel for num_threads(num_of_threads)
			for (j = 0; j <= subseq_count + h; j++)
			{
				BitMap[i][j] = BitMap[i][j] && (bsf[i] > LB[i][j][0]) && (bsf[i] > LB[i][j][3]);
			}
			cnt = 0;

			if(cut_off_neighbors==true)
			{
				for (j = 0; j <= subseq_count + h; j++)
				{
					if (BitMap[i][j])
					{
						if ((cnt > 0) && ((j - CandIndex[cnt-1]) < l_par))
						{
 							if (LB[i][j][3] < LB[i][CandIndex[cnt-1]][3])
							{
 								CandIndex[cnt-1] = j;
							} 
						}
						else
						{
							CandIndex[cnt] = j;
							cnt++;
						}
					}
				}
			}
			else
			{			
				for (j = 0; j <= subseq_count + h; j++)
				{
					if (BitMap[i][j])
					{
						CandIndex[cnt] = j;
						cnt++;
					}
				}
			}
			
			PRF_FINISH(finish5);
			if (cnt == 0) break;

			PRF_START(start6);
			s = 1;
			while(true)
			{
				left = num_of_threads*(s-1);
				right = min(cnt, num_of_threads*s);
				cur_dist = bsf[i];
				
				#pragma omp parallel for num_threads(num_of_threads)
				for (j = left; j < right; j++)
				{
					N[i][CandIndex[j]] = DTW(Q[i], SetC[i][CandIndex[j]], r_par);
					BitMap[i][CandIndex[j]] = false;
				}
				
				for (j = left; j < right; j++)
					if (cur_dist > N[i][CandIndex[j]]) cur_dist = N[i][CandIndex[j]];
								
				s+=1;
				if ((right==cnt)||(bsf[i]>cur_dist)) 
				{
					bsf[i] = min(bsf[i], cur_dist);
					break;
				}
			}
			PRF_FINISH(finish6);
		}
	}
	return 0;
}

void verticalOverlap(int h, bool corr)
{
	int i ,j;
		
	for (j = 0; j < d_par; j++)
	{
		#pragma omp parallel for num_threads(num_of_threads)
		for (i = 0; i <= subseq_count + h; i++)
		{
			if (N[j][i] < inf_val)
			{
				if(corr == true)
				{
					RANK[i] = RANK[i] + ((1 / (N[j][i]+itemtype_epsilon)) * ((d_par - j) / d_par));
				}
				else
				{
					RANK[i] = RANK[i] + (1 / (N[j][i]+itemtype_epsilon));
				}
			}
		}
	}
}


__host__ __device__ itemType LbKim(itemType* Q, itemType* C)
{
	return (dist(Q[0], C[0]) + dist(Q[l_par-1], C[l_par-1]));
}

__host__ __device__ itemType LbKeogh(itemType* Q, itemType* C, int r)
{
	itemType lb_dist = 0;
	itemType min_lokal, max_lokal;

	for (int i = 0; i < l_par; i++)
	{ 
		max_lokal = Q[i];
		min_lokal = Q[i];
		for (int j = max((i - r),0); j <= min((i + r),(l_par-1)); j++)
		{
			max_lokal = max(max_lokal, Q[j]);
			min_lokal = min(min_lokal, Q[j]);
		}
			
		if (C[i] > max_lokal)
		{
			lb_dist = lb_dist + dist(C[i], max_lokal);
		}
		else if (C[i] < min_lokal)
		{
			lb_dist = lb_dist + dist(C[i], min_lokal);
		}
	}
	return lb_dist;
}   

__global__ void LB_calc_kernel(itemType *Q, itemType *SetC, itemType *LB, int N, int idx_max)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int idy = threadIdx.y + blockIdx.y * blockDim.y;
 
 	if(idx<idx_max)
    {

        int s_pos_y=idy*l_par*N;
        int s_pos_x=idx*l_par;
        int start = idy*num_lb*N + idx*num_lb;      
        
        LB[start+0] = LbKim(&Q[idy*l_par],&SetC[s_pos_y + s_pos_x]);
        LB[start+1] = LbKeogh(&Q[idy*l_par],&SetC[s_pos_y + s_pos_x],r_par);
        LB[start+2] = LbKeogh(&SetC[s_pos_y + s_pos_x],&Q[idy*l_par],r_par);
        LB[start+3] = max(LB[start+1], LB[start+2]);
    }
}

