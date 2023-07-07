#include "imputeValue.h"

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

itemType impute(itemType* S, itemType** R, int h)
{

	itemType s = 0;
	PRF_START(start2);
	FillData(S, R, h);
	PRF_FINISH(finish2);
	PRF_START(start3);
	CslcZnorm(h);
	PRF_FINISH(finish3);

	CalcLB(h);
	
	
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

int CalcLB(int h)
{
	int i,j,s,cnt,left,right;
	itemType cur_dist;
	PRF_START(start4);
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