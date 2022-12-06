#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
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

vector<TT> tma(const vector<TT> &h, const vector<TT> &g) {
    vector<TT> alpha, betha;
    alpha.reserve(h.size());
    betha.reserve(h.size());

    alpha.emplace_back(-h[1] / (2 * (h[0] + h[1])));
    betha.emplace_back((3 * (g[1] - g[0])) / (2 * (h[0] + h[1])));

    for (int i = 2; i < h.size(); ++i) {
        alpha.emplace_back(-h[i] / (h[i - 1] * alpha.back() + 2 * (h[i - 1] + h[i])));
        betha.emplace_back(
                (3 * (g[i] - g[i - 1]) - h[i - 1] * betha.back()) / (h[i - 1] * alpha.back() + 2 * (h[i - 1] + h[i])));
    }

    vector<TT> c{0.0};
    c.reserve(h.size());
    for (int i = 1; i < h.size() - 1; ++i) {
        c.emplace_back((3 * (g[i] - g[i - 1]) - h[i - 1] * betha[i]) / (2 * (h[i - 1] + h[i]) + h[i - 1] * alpha[i]));
    }
    c.emplace_back(0.0);
    return c;
}

TT
splineInterpolation(const TT leftBorder, const TT rightBorder, const int countCells, double(*function)(const double),
                    const TT point, bool cellType) {
    vector<pair<TT, TT>> cells = (cellType ? genEvenCells(leftBorder, rightBorder, countCells, function)
                                           : genChebyshovCells(leftBorder, rightBorder, countCells, function));

    if (!cellType) {
        std::reverse(cells.begin(), cells.end());
    }

    //TT resValue = 0.0;
    vector<TT> g;
    g.reserve(countCells);
    vector<TT> h;
    h.reserve(countCells);
    for (int i = 1; i < countCells; ++i) {
        h.emplace_back(cells[i].first - cells[i - 1].first);
        g.emplace_back((cells[i].second - cells[i - 1].second) / h.back());
    }
    vector<TT> c = tma(h, g);

    vector<TT> b;
    vector<TT> d;

    for (int i = 0; i < c.size(); ++i) {
        b.emplace_back(g[i] - ((c[i + 1] + 2 * c[i]) * h[i]) / 3);
        d.emplace_back((c[i + 1] - c[i]) / (3 * h[i]));
        //std::cout << c[i] << " " << d.back() << " " << b.back() << " " << h[i] << std::endl;
    }

    vector<TT> y;
    y.emplace_back(cells[0].second);
    for (int i = 0; i < c.size(); ++i) {
        y.emplace_back(y.back() + b[i] * h[i] + c[i] * pow(h[i], 2) + d[i] * pow(h[i], 3));
    }

    for (double i: y) {
        std::cout << i << "; ";
    }

    std::cout << "\n";

    for (int i = 1; i < cells.size(); ++i) {
//        std::cout << y[i] << "; x = [" << cells[i - 1].first << "; " << cells[i].first << "]; y = ["
//                  << cells[i - 1].second << "; " << cells[i].second << "]\n";
        if (cells[i].first >= point && cells[i - 1].first < point) {
            if (point == cells[i].first)
                return y[i];

            return (y[i - 1] + b[i] * (point - cells[i - 1].first) + c[i] * pow((point - cells[i - 1].first), 2) +
                    d[i] * pow((point - cells[i - 1].first), 3));

        }
    }

    if (rightBorder >= point && cells.back().first < point) {
        return (y.back() + b.back() * (point - cells.back().first) + c.back() * pow((point - cells.back().first), 2) +
                d.back() * pow((point - cells.back().first), 3));
        //return y[i];
    }

    if (cells[0].first > point && leftBorder <= point) {
        return (function(leftBorder) + b[0] * (point - leftBorder) + c[0] * pow((point - leftBorder), 2) +
                d[0] * pow((point - leftBorder), 3));
    }
    //std::cerr << "error\n";
    return y[0];
}


TT
getDiff(const TT leftBorder, const TT rightBorder, const int helpCountCells, const int mainCountCells,
        double(*function)(const double), bool cellType,
        TT(*interpolationFunc)(const TT, const TT, const int, double(*function)(const double), const TT, bool)) {
    vector<pair<TT, TT>> helpCells = (cellType ? genEvenCells(leftBorder, rightBorder, helpCountCells, function)
                                               : genChebyshovCells(leftBorder, rightBorder, helpCountCells,
                                                                   function));
    TT maxValue = 0.0;
    for (const auto &helpCell: helpCells) {
        maxValue = std::max(
                std::abs(function1(helpCell.first) -
                         interpolationFunc(-1.0, 1.0, mainCountCells, function1, helpCell.first, true)),
                maxValue);
    }
    return maxValue;
}

void getLagrange() {
    //std::cout << "Lagrange, even cells: " << lagrangeInterpolation(-1.0, 1.0, 3, function1, 0, true);
    //std::cout << "\nLagrange, Chebyshov: " << lagrangeInterpolation(-1.0, 1.0, 3, function1, 0, false);
    //std::cout << "\nnLagrange DIFF: " << getDiff(-1.0, 1.0, 512, 100, function1, false, lagrangeInterpolation);
}

void getSpline() {
    std::cout << "Spline, even cells: " << splineInterpolation(-1.0, 1.0, 3, function1, 1, true);
    std::cout << "\nSpline, Chebyshov: " << splineInterpolation(-1.0, 1.0, 3, function1, -1, false);
    //std::cout << "\nSpline DIFF: " << getDiff(-1.0, 1.0, 512, 100, function1, false, splineInterpolation);
}
