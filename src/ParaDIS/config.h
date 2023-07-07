#ifndef CONFIG_H
#define CONFIG_H
#include <math.h>

#define L_par (22500) //start point
#define lengthsrt (5000)
#define numimputevalue (2500) // count of imputation points

#define cum_error (true) 

//init bsf: 
//0 - random, 
//1 - prev point
//2 - min(max(LB))
#define init_bsf (2) 
#define cut_off_neighbors (true) 

#define d_par (3) 
#define l_par (50) 
#define k_par (3) 

#define r_par (int(floor(l_par/4))) 
#define inf_val (1e+10)
#define itemtype_epsilon (2.2204460492503131e-16)
#define num_lb (4)
#define use_corr (false)//(false)

#define dist(x,y) (((x)-(y))*((x)-(y)))
#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))  

#define subseq_count ((L_par) - (l_par))

typedef float itemType;

#endif // CONFIG_H
