/*+++++++++++++++++++++++++++++++++

Project: ParaDIS (Parallel Algorithm for Imputation of Missing Values in Streaming Time Series)

Source file: ReadData.cpp

Purpose: Read dataset from file

Author(s): Andrey Poluyanov (andrey.poluyanov@gmail.com) and Mikhail Zymbler (mzym@susu.ru)

+++++++++++++++++++++++++++++++++*/

#include "ReadData.h"

void read_timeseries(const char* filename, itemType *S, itemType **R)
{	
	FILE* fileS;
		fileS = fopen(filename, "r"); 
	
	assert(fileS != NULL); 
	
	char line[lengthsrt];

	int i = 0;
	while (fgets(line, lengthsrt, fileS))
	{
		char* tmp = strdup(line);
		
		S[i] = atof(getelem(tmp, 1));
		free(tmp);
		for (int j = 0; j < d_par; j++)
		{
			char* tmp = strdup(line);
			R[j][i] = atof(getelem(tmp, j + 2));
			free(tmp);
		}
		i++;
		if (i >= L_par + numimputevalue)
			break;
		
	}

}


const char* getelem(char* line, int num)
{
	const char* tok;
	for (tok = strtok(line, ";"); tok && *tok; tok = strtok(NULL, ";\n"))
	{
		if (!--num)
		{
			return tok;
		}
	}
	return NULL;
}


