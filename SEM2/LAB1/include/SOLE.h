#ifndef INC_LAB_SOLE_H
#define INC_LAB_SOLE_H

#include <vector>
#include "shared.h"

// const std::pair<TT, TT> range = {0, 3.5};

// test 1
// const std::vector<TT> initPoints = {0.0, -0.991};
// const std::vector<TT> initPoints = {0.0, 0.991};

// test 2
// const std::vector<TT> initPoints = {0.0, 0.991};
// const std::vector<TT> initPoints = {0.0, -0.991};

// test 3
// const std::vector<TT> initPoints = {1.0, 2.0, 3.0};

const std::vector<TT> initPoints = {1.0, 0.0};
const int numOfPoints = 1000;
const std::pair<TT, TT> range = {0.0, 10.0};
const TT step = TT(range.second - range.first) / TT(numOfPoints);

#endif //INC_LAB_SOLE_H
