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
