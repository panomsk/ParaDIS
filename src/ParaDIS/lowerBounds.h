/*+++++++++++++++++++++++++++++++++

Project: ParaDIS (Parallel Algorithm for Imputation of Missing Values in Streaming Time Series)

Source file: lowerBounds.h

Purpose: Lower bounds calculation for two subsequences

Author(s): Andrey Poluyanov (andrey.poluyanov@gmail.com) and Mikhail Zymbler (mzym@susu.ru)

+++++++++++++++++++++++++++++++++*/

#ifndef LOWERBOUNDS_H
#define LOWERBOUNDS_H

#include "config.h"

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>

itemType LbKim(itemType* Q, itemType* C);

itemType LbKeogh(itemType* Q, itemType* C, int r);

#endif // LOWERBOUNDS_H
