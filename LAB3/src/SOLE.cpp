#include <iostream>
#include <string>
#include <cmath>
#include "shared.h"

using std::vector;
using std::string;
using std::pair;

TT function1(const TT x) {
    return pow(x, 2);
}

vector<pair<TT, TT>>
genEvenCells(const TT leftBorder, const TT rightBorder, const int countCells, double(*function)(const double)) {
    TT step = (rightBorder - leftBorder) / (countCells - 1);
    vector<pair<TT, TT>> cell;
    cell.reserve(countCells);
    for (int i = 0; i < countCells; ++i) {
//        cell[i].first = leftBorder + step * i;
//        cell[i].second = function(cell[i].first);
        cell.emplace_back(leftBorder + step * i, function(leftBorder + step * i));
        //std::cout << "point " << cell[i].first << "; value " << cell[i].second << "\n";
    }
    return cell;
}

vector<pair<TT, TT>>
genChebyshovCells(const TT leftBorder, const TT rightBorder, const int countCells, double(*function)(const double)) {
    TT multiplayer = (rightBorder - leftBorder) / 2;
    TT addition = (rightBorder + leftBorder) / 2;
    vector<pair<TT, TT>> cell;
    cell.reserve(countCells);
    for (int i = 0; i < countCells; ++i) {
        TT temp = addition + multiplayer * cos((2 * i + 1) * M_PI / (2 * (countCells + 1)));
        cell.emplace_back(temp, function(temp));
        //std::cout << "point " << cell[i].first << "; value " << cell[i].second << "\n";
    }
    return cell;
}

TT calcC(const vector<pair<TT, TT>> &cells, const TT point, const int splitPoint) {
    TT c = 1.0;
    //std::cout << "\npoint: " << point << "; split: " << splitPoint << "\n";
    for (int i = 0; i < splitPoint; ++i) {
        c *= (point - cells[i].first) / (cells[splitPoint].first - cells[i].first);
        //std::cout << "<- C: " << c << "; ";
    }
    for (int i = splitPoint + 1; i < cells.size(); ++i) {
        c *= (point - cells[i].first) / (cells[splitPoint].first - cells[i].first);
        //std::cout << "C B: " << c;
    }
    return c;
}

//TT calc diff(const vector<pair<TT, TT>> &cells, double(*function)(const double)){
//    for (int i = 0; i <)
//}

TT
lagrangeInterpolation(const TT leftBorder, const TT rightBorder, const int countCells, double(*function)(const double),
                      const TT point, bool cellType) {
    vector<pair<TT, TT>> cells = (cellType ? genEvenCells(leftBorder, rightBorder, countCells, function)
                                           : genChebyshovCells(leftBorder, rightBorder, countCells, function));
    TT resValue = 0.0;
    for (int i = 0; i < countCells; ++i) {
        resValue += calcC(cells, point, i) * cells[i].second;
        //std::cout << "cVal " << calcC(cells, point, i) << "; value " << cells[i].second << "\n";
    }
    return resValue;
}

TT
getDiff(const TT leftBorder, const TT rightBorder, const int helpCountCells, const int mainCountCells,
        double(*function)(const double), bool cellType) {
    vector<pair<TT, TT>> helpCells = (cellType ? genEvenCells(leftBorder, rightBorder, helpCountCells, function)
                                               : genChebyshovCells(leftBorder, rightBorder, helpCountCells,
                                                                   function));
    TT maxValue = 0.0;
    for (const auto &helpCell: helpCells) {
        maxValue = std::max(
                std::abs(function1(helpCell.first) -
                         lagrangeInterpolation(-1.0, 1.0, mainCountCells, function1, helpCell.first, true)),
                maxValue);

//        TT temp = lagrangeInterpolation(-1.0, 1.0, mainCountCells, function1, helpCell.first, true);
//        for (const auto &secondHelp: helpCells) {
//            maxValue = std::max(
//                    std::abs(temp -
//                             lagrangeInterpolation(-1.0, 1.0, mainCountCells, function1, secondHelp.first, true)),
//                    maxValue);
//        }
    }
    return maxValue;
}

void getLagrange() {
    //std::cout << "Lagrange, even cells: " << lagrangeInterpolation(-1.0, 1.0, 3, function1, 0, true);
    //std::cout << "\nLagrange, Chebyshov: " << lagrangeInterpolation(-1.0, 1.0, 3, function1, 0, false);
    std::cout << "\nDIFF: " << getDiff(-1.0, 1.0, 512, 100, function1, false);
}
