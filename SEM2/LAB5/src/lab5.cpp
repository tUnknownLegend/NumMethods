#include <iostream>
#include <fstream>
#include "test.h"
#include "SOLE.h"

using namespace std;

void solve_quadrature(const std::string &path, const size_t _sections) {
    TT h = (area.back() - area.front()) / _sections;
    size_t number_of_points = _sections + 1;

    std::vector<TT> grid;
    grid.reserve(number_of_points);
    std::vector<TT> right_part;
    right_part.reserve(number_of_points);
    std::vector<TT> coef_a(number_of_points, h);
    coef_a.front() = h / 2;
    coef_a.back() = h / 2;

    for (size_t i = 0; i < number_of_points; ++i) {
        grid.push_back(area.front() + i * h);
        right_part.push_back(f(grid[i]));
    };

    vector<vector<TT>> matrix(number_of_points, vector<TT>(number_of_points));

    for (size_t i = 0; i < number_of_points; ++i) {
        for (size_t k = 0; k < number_of_points; ++k)
            matrix[i][k] = -lambda * K(grid[i], grid[k]) * coef_a[k];
        matrix[i][i] += 1;
    }

    const auto res = CalcGaussMethod(matrix, right_part);

    std::ofstream file;
    file.open(path);
    for (size_t i = 0; i < res.size(); ++i) file << grid[i] << " " << res[i] << "\n";
    file.close();
}


void solve_simple_iterations(const std::string &path, const size_t _sections, const double eps,
                             int iterations) {
    std::ofstream file;
    file.open(path);
    TT h = (area.back() - area.front()) / _sections;
    size_t number_of_points = _sections + 1;

    std::ofstream errors;
    errors.open("../output/errors.dat");
    std::vector<TT> result;
    result.reserve(number_of_points);
    for (size_t i = 0; i < number_of_points; ++i) result.push_back(f(area.front() + i * h));
    std::vector<TT> iterated(number_of_points, TT(0));
    const int _iterations = iterations;
    while (iterations > 0) {
        for (size_t i = 0; i < number_of_points; ++i) {
            const TT x_i = area.front() + i * h;
            const TT f_i = f(x_i);
            const TT c_i = [&h, &number_of_points, &x_i, &result]() -> TT {
                TT sum = TT(0);
                for (size_t k = 0; k < number_of_points; ++k) {
                    const TT x_k = area.front() + k * h;
                    const TT a_k = (k == 0 || k == number_of_points - 1) ? h / 2 : h;
                    sum += lambda * K(x_i, x_k) * result[k] * a_k;
                }
                return sum;
            }();
            iterated[i] = c_i + f_i;
        }
        const TT _norm = norm1Vector(vectorOperation(result, iterated, '-'));
        errors << _iterations - iterations + 1 << " " << _norm << "\n";
        result = iterated;
        --iterations;
        for (size_t i = 0; i < result.size(); ++i) {
            file << area.front() + i * h << " " << result[i] << " ";
        }
        file << "\n";
    }
    errors.close();

    file.close();
}

void solve_degenerate(const std::string &path, const size_t sections) {
    TT h = (area.back() - area.front()) / sections;
    size_t size = phi.size();
    vector<vector<TT>> matrix(size, vector<TT>(size));
    std::vector<TT> right_part;
    right_part.reserve(size);
    std::vector<TT> x;
    x.reserve(sections);
    for (size_t i = 0; i < sections + 1; ++i) x.push_back(area.front() + (i) * h);

    auto integrate = [sections, h, x](std::function<TT(TT)> f, std::function<TT(TT)> g) -> TT {
        TT sum = 0;
        for (size_t i = 0; i < sections; ++i)
            sum += (f(x[i]) * g(x[i]) + f(x[i + 1]) * g(x[i + 1])) / 2;
        return (sum * h);
    };

    for (size_t i = 0; i < size; ++i) {
        right_part.push_back(integrate(psi[i], f));
        for (size_t j = 0; j < size; ++j) {
            matrix[i][j] = -lambda * integrate(psi[i], phi[j]);
        }
        matrix[i][i] += 1;
    }

    const auto result = CalcGaussMethod(matrix, right_part);

    std::ofstream file;
    file.open(path);
    for (const auto &i: result) {
        file << i << "\n";
    }
    file.close();

}


TT solve_singular(const std::string &path, const size_t sections) {
    const TT length = 2 * M_PI / sections;
    std::vector<std::array<TT, 2>> k;
    k.reserve(sections);
    std::vector<std::array<TT, 2>> c;
    c.reserve(sections);
    std::vector<std::array<TT, 2>> n;
    n.reserve(sections);
    std::vector<TT> right_part;
    right_part.reserve(sections + 1);
    vector<vector<TT>> matrix(sections + 1, vector<TT>(sections + 1, 1.0));

    for (size_t i = 1; i < sections + 1; ++i) {
        const TT var1 = length * (i - 0.5);
        const TT var2 = length * (i - 1);

        k.push_back({cos(var1), sin(var1)});
        c.push_back({cos(var2), sin(var2)});
        n.push_back({cos(var1), sin(var1)});
        right_part.push_back(F(var1));
    }
    right_part.push_back(TT{0});

    const auto scalar = [](const std::array<TT, 2> &a, const std::array<TT, 2> &b) -> TT {
        return (a[0] * b[0] + a[1] * b[1]);
    };

    for (size_t i = 0; i < sections; ++i) {
        for (size_t j = 0; j < sections; ++j) {
            matrix[i][j] *= scalar(Q(k[i], c[j]), n[i]);
        }
        matrix[i][sections] = TT{1};
    }
    matrix[sections][sections] = TT{0};

    const auto result = CalcGaussMethod(matrix, right_part);

    std::ofstream file;
    file.open(path);
    for (size_t i = 0; i < result.size(); ++i) {
        file << c[i][0] << " " << c[i][1] << " " << result[i] << "\n";
    }
    file.close();
    return result.back();
}
