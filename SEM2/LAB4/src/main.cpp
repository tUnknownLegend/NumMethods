#include "lab4.h"
#include <cmath>
#include <vector>
#include <utility>
#include "shared.h"
//#include "SOLE.h"


void depGraphData() {
    const unsigned int numOfPoints = 100000;

    std::array<size_t, 2> N = {10, 10};
    const std::string path = "../output/temp.dat";

    std::vector<std::pair<TT, TT>> graphData;
    graphData.reserve(numOfPoints);

    for (unsigned int i = 1; i < numOfPoints / 100; ++i) {
        graphData.emplace_back((TT) i / (TT) numOfPoints, solve(path, N,
                                                                (TT) i / (TT) numOfPoints, 1e-8));
    }

    outputPair(graphData, "../output/graphData.txt");
}


int main() {
    const std::string path = "../output/solution.dat";

    // Количество разбиений (x1 / x2 / время)
    std::array<size_t, 2> N = {100, 100};

    solve(path, N, 0.001, 1e-8);

    //  depGraphData();

    return 0;
}
