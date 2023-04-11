#ifndef INC_LAB_SOLE_H
#define INC_LAB_SOLE_H

#include <cmath>
#include "shared.h"

const unsigned int lh = 1;
const unsigned int lv = 1;

/* steps */
// horizontal step
const TT h = 0.1;
// vertical step
const TT tao = 0.01;
/**/
const int n = round(lh / h) + 1;
const int k = round(lv / tao) + 1;

/**/

const TT a = 1.0;


#endif //INC_LAB_SOLE_H
