#ifndef INC_LAB1_QR_H
#define INC_LAB1_QR_H

#include <vector>
#include "shared.h"
//#define TT double

std::vector<TT> CalcQRmethod(std::vector<std::vector<TT>>& matrix, std::vector<TT> rightVect);

std::vector<TT> getQR();

#endif //INC_LAB1_QR_H
