#ifndef INC_LAB2_SOLE_H
#define INC_LAB2_SOLE_H

#include <vector>
#include "shared.h"

const std::vector<TT> initPoints = {1.0, 0.0};
const int numOfPoints = 10;
const std::pair<TT, TT> range = {0, 1};
const TT step = TT(range.second - range.first) / TT(numOfPoints) ;

#endif //INC_LAB2_SOLE_H
