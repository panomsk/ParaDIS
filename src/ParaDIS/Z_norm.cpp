/*+++++++++++++++++++++++++++++++++

Project: ParaDIS (Parallel Algorithm for Imputation of Missing Values in Streaming Time Series)

Source file: Z_norm.cpp

Purpose: Z-normalization of subsequence

Author(s): Andrey Poluyanov (andrey.poluyanov@gmail.com) and Mikhail Zymbler (mzym@susu.ru)

+++++++++++++++++++++++++++++++++*/

#include "Z_norm.h"

void Z_norm(itemType* query, int lenQuery)
{
	itemType mean = 0;
	for (int i = 0; i < lenQuery; i++)
		mean = mean + query[i];
	mean = mean / lenQuery;

	itemType variance = 0;
	for (int i = 0; i < lenQuery; i++)
		variance = variance + (query[i] - mean)*(query[i] - mean);
	variance = variance / lenQuery;
	
	itemType sdiviation = sqrt(variance);
	for (int i = 0; i < lenQuery; i++)
		query[i] = (query[i] - mean) / sdiviation;
}


