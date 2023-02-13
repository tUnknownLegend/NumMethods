#include <iostream>
#include <vector>
#include "shared.h"

using namespace std;

template<typename T>
using matrix = vector<vector<T>>;

template<typename T, typename F>
auto derivative(F& fun, T eps) {
    return [&fun, eps](veT x) {
        return (fun(x + eps) - fun(x)) / eps;
    };
}

vector<TT> Newton(function<TT(vector<TT> &)> &fun, TT diff, pair<TT, TT> range,
          const vector<TT> &initPoints, TT eps, size_t &iterCount) {
    vector<TT> xk = initPoints;
    TT diffVal;
    do {
        vector<TT> newXk = xk - fun(xk) / diff(xk);
        diffVal = normInfVector(vectorOperation(newXk, xk, '-'));
        xk = newXk;
        ++iterCount;
    } while (diffVal > eps);
    return xk;
}
