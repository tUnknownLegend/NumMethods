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

    std::array<size_t, 2> N = {10, 10}; //Количество разбиений (x1 / x2 / время)

    solve(path, N, 0.001, 1e-8);

    depGraphData();

//    const std::string solution2x = "../output/solution 2x.dat";
//    solve(solution2x, {20, 20}, 0.0005, 1e-4);
//
//    const std::string solution_1x1y = "../output/solution_1x1y.dat";
//    const std::string solution_2x1y = "../output/solution_2x1y.dat";
//    const std::string solution_3x1y = "../output/solution_3x1y.dat";
//    const std::string solution_1x2y = "../output/solution_1x2y.dat";
//    const std::string solution_2x2y = "../output/solution_2x2y.dat";
//    const std::string solution_3x2y = "../output/solution_3x2y.dat";
//    const std::string solution_1x3y = "../output/solution_1x3y.dat";
//    const std::string solution_2x3y = "../output/solution_2x3y.dat";
//    const std::string solution_3x3y = "../output/solution_3x3y.dat";
//
//    solve(solution_1x1y, {5, 5}, 0.001, 1e-4);
//    solve(solution_2x1y, {10, 5}, 0.001, 1e-4);
//    solve(solution_3x1y, {20, 5}, 0.001, 1e-4);
//    solve(solution_1x2y, {5, 10}, 0.001, 1e-4);
//    solve(solution_2x2y, {10, 10}, 0.001, 1e-4);
//    solve(solution_3x2y, {20, 10}, 0.001, 1e-4);
//    solve(solution_1x3y, {5, 20}, 0.001, 1e-4);
//    solve(solution_2x3y, {10, 20}, 0.001, 1e-4);
//    solve(solution_3x3y, {20, 20}, 0.001, 1e-4);

    return 0;
}
