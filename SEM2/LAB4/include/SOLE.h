#ifndef INC_LAB_SOLE_H
#define INC_LAB_SOLE_H

#include <cmath>
#include "shared.h"

const unsigned int lh = 4;
const unsigned int lv = 4;

const unsigned int leftBorder = -2;

///* steps */
//// horizontal step
//const TT h = 0.01;
//// vertical step
//const TT tao = 0.01;
///**/
//const int n = round(lh / h) + 1;
//const int k = round(lv / tao) + 1;
//
//const TT a = 1.0;

std::vector<TT> tma(std::vector<TT> vecA, std::vector<TT> vecB, std::vector<TT> vecC, std::vector<TT> vecD);

#endif //INC_LAB_SOLE_H
