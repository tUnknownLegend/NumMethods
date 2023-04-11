#ifndef INC_LAB_SOLE_H
#define INC_LAB_SOLE_H

#include <cmath>
#include "shared.h"

/* equation params */
// linear density
const TT ro = 1.0;
// specific heat capacity
const TT c = 1.0;

/* stick data */
// pivot temperature point of stick
const TT t0 = 1.0;
// length
const TT l = 1.0;

const TT Q = 10;

/* steps */
// horizontal step
const TT h = 0.2;
//const TT h = 0.08;
// vertical step
const TT tao = 0.2;
//const TT tao = 0.002;
/**/
const int N = round(l / h - 1);
const int k = round(t0 / tao);

const TT sigma = 0.0;

// K calc
const TT k1 = 0.1;
const TT k2 = 1.5;
const TT x1 = 0.25;
const TT x2 = 0.5;

#endif //INC_LAB_SOLE_H
