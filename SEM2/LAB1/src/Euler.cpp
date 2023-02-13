#include <cmath>
#include <algorithm>
#include <utility>
#include "shared.h"
#include "nonLinearSolve.h"

using namespace  std;

const vector<TT> initCells = {0.0, 1.0};
const TT step = 0.1;
const pair<TT, TT> range = {0, 1};

TT calcFunc(const vector<TT>& args) {
    return pow(args.at(0), 2) - 2 * args.at(1);
}

vector<TT> calcSingleCellSimple(const vector<TT>& prevPoints) {
    vector<TT> currPoints;
    currPoints.reserve(prevPoints.size());

    currPoints.push_back(prevPoints.at(0) + step);
    for (size_t i = 1; i < prevPoints.size() - 1; ++i) {
        currPoints.push_back(prevPoints.at(1) + step * calcFunc(prevPoints));
    }
    currPoints.push_back(calcFunc(currPoints));
    return currPoints;
}

vector<TT> calcSingleCellExplicit(const vector<TT>& prevPoints) {
    vector<TT> currPoints;
    currPoints.reserve(prevPoints.size());



    currPoints.push_back(prevPoints.at(0) + step);
    for (size_t i = 1; i < prevPoints.size() - 1; ++i) {
        currPoints.push_back(prevPoints.at(1) + step * calcFunc(prevPoints));
    }
    currPoints.push_back(calcFunc(currPoints));
    return currPoints;
}

vector<vector<TT>> calcCell(const function<vector<TT>(const vector<TT>&)>& singleCell) {
    int const cellSize = ceil((range.second - range.first) / step);
    vector<vector<TT>> cells;
    cells.reserve(cellSize);

    vector<TT> currCells( initCells );
    currCells.push_back(calcFunc(initCells));

    for(int i = 0; i < cellSize; ++i) {
        cells[i].reserve(currCells.size());
        cells.push_back(singleCell(currCells));
        currCells = cells.back();
    }
    return cells;
}

void Euler() {
    outputOnTheScreenMatrix(calcCell(calcSingleCellSimple));
}

void ExplicitEuler() {
    outputOnTheScreenMatrix(calcCell(calcSingleCellExplicit));
}

