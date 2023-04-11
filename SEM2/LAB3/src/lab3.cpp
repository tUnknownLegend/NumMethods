#include "SOLE.h"
#include <iostream>

using namespace std;

TT f(const TT x) {
    // test 1
    // return sin(M_PI * x);

    // test 2
//    return x - pow(x, 2);

    // test 13
    return (pow(x, 2) + 1) / 2;
}

TT Fxx(const TT x) {
    // test 1
    // return -M_PI * M_PI * sin(M_PI * x);

    // test 2
//    return -2.0;

    // test 13
    return 1.0;
}

TT g(const TT x) {
//    return 0.0;

    // test 13
    return x * sin(2 * x);
}

TT phi(const TT t) {
//    return 0.0;

    // test 13
//    return 0.5 + 3 * t;
}

TT rho(const TT t) {
//    return 0.0;

    // test 13
//    return 1.0;
}

TT preciseSolution(const int test,
                   const unsigned int i, const unsigned int j) {
    TT solution = 0.0;
    const TT eps = 1e-6;
    const TT K = (sqrt(2 / (pow(M_PI, 2) * eps)) - 1) / 2;

    switch (test) {
        case 1:
            solution = sin(M_PI * i * h) * cos(M_PI * j * tao);
            break;
        case 2:
            for (unsigned int l = 0; l < K; ++l) {
                solution += 1 / pow((2 * l + 1), 3) *
                            sin(2 * l + 1) * M_PI + i * h *
                                                    cos(2 * l + 1) * M_PI * j * tao;
            }
            solution *= 8 / pow(M_PI, 3);
            break;
        default:
            cerr << "no such test: " << test;
    }

    return solution;
}

vector<vector<TT>> approxInitialCondition() {
    vector<vector<TT>> result(k, vector<TT>(n, 0.0));

    for (unsigned int i = 0; i < n; ++i) {
        result[0][i] = f(i * h);
    }

    for (unsigned int j = 0; j < k; ++j) {
        result[j][0] = phi(j * tao);
        result[j][n - 1] = rho(j * tao);
    }

    return result;
}

void getFirstStep(vector<vector<TT>> &result) {
    for (unsigned int i = 0; i < n; ++i) {
        result[1][i] = result[0][i] + tao * g(i * h) +
                       pow(a, 2) * pow(tao, 2) / 2 *
                       Fxx(i * h);
    }
}

void getFurtherSteps(vector<vector<TT>> &result) {
    for (unsigned int j = 1; j < k - 1; ++j) {
        for (unsigned int i = 1; i < n - 1; ++i) {
            result[j + 1][i] = pow(a, 2) / pow(h, 2) *
                               (result[j][i + 1] - 2 * result[j][i] + result[j][i - 1]) *
                               pow(tao, 2) -
                               (-2 * result[j][i] + result[j - 1][i]);
        }
    }
}

vector<vector<TT>> crossScheme() {
    vector<vector<TT>> result = approxInitialCondition();

    getFirstStep(result);
    getFurtherSteps(result);

    return result;
}

TT getDiff(const vector<vector<TT>> &result) {
    TT diff = result[0][0];

    for (unsigned int j = 0; j < k; ++j) {
        for (unsigned int i = 0; i < n; ++i) {
            diff = max(diff, abs(
                    preciseSolution(2, i, j) -
                    result[j][i]));
        }
    }

    return diff;
}

void getLab3Result() {
    vector<vector<TT>> result = crossScheme();

    outputMatrix(result);
    cout << "diff: " << getDiff(result) << "\n";
}
