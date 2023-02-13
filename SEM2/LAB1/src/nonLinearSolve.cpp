#include <iostream>
#include <vector>
#include "shared.h"

using namespace std;

template<typename T>
using matrix = vector<vector<T>>;

auto derivative(function<TT(vector<TT> &)> &fun, TT eps) {
    return [&fun, eps](vector<TT> &x) {
        vector<TT> xDiff(x);
        vectorDigit(eps, xDiff, '+');
        return ((fun(xDiff) - fun(x)) / eps);
    };
}

vector<TT> Newton(function<TT(vector<TT> &)> &fun,
                  const vector<TT> &initPoints, TT eps, size_t &iterCount) {
    vector<TT> currPoints = initPoints;
    auto diff = derivative(fun, COMPARE_RATE);
    TT diffVal;
    do {
        TT der = diff(currPoints);
        vectorDigit(, der, '/');
        vector<TT> newXk = vectorOperation(fun(currPoints) / der - der, der, '-');
        diffVal = normInfVector(vectorOperation(newXk, currPoints, '-'));
        currPoints = newXk;
        ++iterCount;
    } while (diffVal > eps);
    return currPoints;
}
