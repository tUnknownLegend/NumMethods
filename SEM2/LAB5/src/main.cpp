#include "lab5.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include "shared.h"

std::vector<std::pair<TT, TT>> getRstats(const std::string &path, const std::string &format) {
    std::vector<std::pair<TT, TT>> Rstat;

    for (unsigned int i = 10; i <= 1000; i += 10) {
        Rstat.push_back({i, solve_singular(path + "singular" + format, i)});
    }
    return Rstat;
}


int main() {
    const std::string path = "../output/";
    const std::string format = ".dat";

    //Количество разбиений
    const size_t N = 100;

//    solve_quadrature(path + "quad" + format, N);

//    solve_simple_iterations(path + "simple" + format, N, 1e-10, 100);

//    solve_degenerate(path + "degenerate" + format, 80);

//    solve_singular(path + "singular" + format, N);

    outputPair(getRstats(path, format), path + "Rstats" + format);

    return 0;
}
