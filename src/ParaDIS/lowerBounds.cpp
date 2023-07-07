#include "lowerBounds.h"


itemType LbKim(itemType* Q, itemType* C)
{
	return (dist(Q[0], C[0]) + dist(Q[l_par-1], C[l_par-1]));
}

itemType LbKeogh(itemType* Q, itemType* C, int r)
{
	itemType lb_dist = 0;
	itemType upper[l_par], lower[l_par];

	for (int i = 0; i < l_par; i++)
	{ 
		itemType max_lokal = Q[i];
		itemType min_lokal = Q[i];
		//!! исправлено значение верхней границы 
		for (int j = max((i - r),0); j <= min((i + r),(l_par-1)); j++)
		{
			max_lokal = max(max_lokal, Q[j]);
			min_lokal = min(min_lokal, Q[j]);
		}
		//!! вынесено во внешний цикл
		upper[i] = max_lokal;
		lower[i] = min_lokal;
	}
		
	for (int i = 0; i < l_par; i++)
	{
		if (C[i] > upper[i])
		{
			lb_dist = lb_dist + dist(C[i], upper[i]);
		}
		else if (C[i] < lower[i])
		{
			lb_dist = lb_dist + dist(C[i], lower[i]);
		}
	}
	//printf_s("%lf\n", lb_dist);
	return lb_dist;
}   