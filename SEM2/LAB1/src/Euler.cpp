#include "shared.h"
#include "nonLinearSolve.h"
#include "SOLE.h"

using namespace std;

vector<vector<TT>> explicitEuler(const vector<TT> &cond, int n) {
    vector<vector<TT>> res(n, vector<TT>(cond.size()));
    res[0] = cond;

    for (int i = 0; i < n - 1; ++i) {
        res[i + 1] = vectorOperation(res[i],
                                     vectorDigit(step, f(res[i]), '*'), '+');
    }
    return res;
}

vector<vector<TT>> implicitEuler(const vector<TT>& cond, const int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    for (int i = 0; i < n - 1; i++) {
        y[i + 1] = Newton("Euler", cond.size(), y[i]);
    }
    //cout << "Решение в k-м узле: " << y[t][0] << ", " << y[t][1] << endl;
    return y;
}

// Симметричная схема
//vector<vector<TT>> symmetric(vector<TT> cond, TT h, int n, int dim) {
//    std::string flag = "Symmetric";
//    vector<vector<TT>> y(n, vector<TT>(dim));
//    vector<TT> temp(dim);
//    y[0] = cond;
//    //int t = 1200;
//
//    for (int i = 0; i < n - 1; i++)
//    {
//        temp = y[i] + h * f(i * h, y[i], dim);
//        y[i + 1] = Newton(flag, dim, temp, y[i], y[i + 1], h, (i + 1) * h);
//    }
//    //cout << "Решение в k-м узле: " << y[t][0] << ", " << y[t][1] << endl;
//    return y;
//}

//TT calcFunc(const vector<TT> &args) {
//    return pow(args.at(0), 2) - 2 * args.at(1);
//}
//
//vector<TT> calcSingleCellSimple(const vector<TT> &prevPoints) {
//    vector<TT> currPoints;
//    currPoints.reserve(prevPoints.size());
//
//    currPoints.push_back(prevPoints.at(0) + step);
//    for (size_t i = 1; i < prevPoints.size() - 1; ++i) {
//        currPoints.push_back(prevPoints.at(i) + step * calcFunc(prevPoints));
//    }
//    currPoints.push_back(calcFunc(currPoints));
//    return currPoints;
//}
//
//vector<TT> calcSingleCellImplicit(const vector<TT> &prevPoints) {
//    vector<TT> currPoints;
//    currPoints.reserve(prevPoints.size());
//
//    //nonLinearSolve(calcFunc, initPoints, COMPARE_RATE, {});
//
//    currPoints.push_back(prevPoints.at(0) + step);
//    for (size_t i = 1; i < prevPoints.size() - 1; ++i) {
//        currPoints.push_back(prevPoints.at(1) + step * calcFunc(prevPoints));
//    }
//    currPoints.push_back(calcFunc(currPoints));
//    return currPoints;
//}
//
//vector<vector<TT>> calcCell(const function<vector<TT>(const vector<TT> &)> &singleCell) {
//    int const cellSize = ceil((range.second - range.first) / step);
//    vector<vector<TT>> cells;
//    cells.reserve(cellSize);
//
//    vector<TT> currCells(initPoints);
//    currCells.push_back(calcFunc(initPoints));
//
//    for (int i = 0; i < cellSize; ++i) {
//        cells[i].reserve(currCells.size());
//        cells.push_back(singleCell(currCells));
//        currCells = cells.back();
//    }
//    return cells;
//}

void ImplicitEuler() {
    outputOnTheScreenMatrix(explicitEuler(initPoints, sizeN));
}

void ExplicitEuler() {
    outputOnTheScreenMatrix(explicitEuler(initPoints, sizeN));
}

