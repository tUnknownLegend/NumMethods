#include "SOLE.h"
#include <iostream>

using namespace std;

TT f(const TT x) {
    return sin(M_PI * x);
}

TT Fxx(const TT x) {
    return -M_PI * M_PI * sin(M_PI * x);
}

TT g(const TT x) {
    return 0.0;
}

TT phi(const TT t) {
    return 0.0;
}

vector<vector<TT>> approxInitialCondition() {
    vector<vector<TT>> result(k, vector<TT>(n, 0.0));

    for (unsigned int i = 0; i < n; ++i) {
        result[0][i] = f(i * h);
    }

    for (unsigned int j = 0; j < k; ++j) {
        result[j][0] = phi(j * tao);
        result[j][n - 1] = phi(j * tao);
    }

    return result;
}

void getFirstStep(vector<vector<TT>> &result) {
    for (unsigned int i = 1; i < n - 1; ++i) {
        result[1][i] = result[0][i] + tao * g(i * h) +
                       pow(a, 2) * pow(tao, 2) / 2 *
                       Fxx(i * h);
    }
}

void getFurtherSteps(vector<vector<TT>> &result) {
    for (unsigned int j = 2; j < k; ++j) {
        for (unsigned int i = 1; i < n - 1; ++i) {
            result[j][i] = pow(a, 2) / pow(h, 2) *
                           (result[j - 1][i + 1] - 2 * result[j - 1][i] + result[j - 1][i - 1]) *
                           pow(tao, 2) -
                           (-2 * result[j - 1][i] + result[j - 2][i]);
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
                    sin(M_PI * i * h) * cos(M_PI * j * tao) -
                    result[j][i]));
        }
    }

    return diff;
}

void getLab3Result() {
    vector<vector<TT>> result = crossScheme();

    outputMatrix(result);
//    outputOnTheScreenMatrix(result);

    cout << "diff: " << getDiff(result) << "\n";
}
