#ifndef LAB3_LAB3_H
#define LAB3_LAB3_H

#include "shared.h"

void solve_quadrature(const std::string &path, size_t _sections);

void solve_simple_iterations(const std::string &path, size_t _sections,
                             double eps, int iterations);

void solve_degenerate(const std::string &path, size_t sections);

TT solve_singular(const std::string &path, size_t sections);

#endif //LAB3_LAB3_H
