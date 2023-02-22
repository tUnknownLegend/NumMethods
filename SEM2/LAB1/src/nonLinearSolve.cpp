#include <vector>
#include <string>
#include "shared.h"
#include "SOLE.h"

using namespace std;

vector<TT> f(const vector<TT> &x) {
    // example 1
//    return {2 * x[0] + x[1] * x[1] - 1, 6 * x[0] - x[1] * x[1] + 1};
    // example 2
//    return {1 - x[0] * x[0] - x[1] * x[1], 2 * x[0]};
    // example 3
//    return {10 * (x[1] - x[0]), x[0] * (28 - x[2]) - x[1], x[0] * x[1] - 8 / 3 * x[2]};
    // pendulum
    return {2 * x[1], -2 / 0.3 * x[0]};
}

vector<TT> calcMethod(const size_t dim, const string &method, const vector<TT> &x0, const vector<TT> &y0) {
    vector<TT> ans(dim);

    if (method == "Euler") {
        for (size_t i = 0; i < dim; i++) {
            ans[i] = x0[i] - y0[i] - step * f(x0)[i];
        }
    }

    if (method == "Symmetric")
        for (size_t i = 0; i < dim; i++) {
            ans[i] = x0[i] - y0[i] - step / 2 * (f(x0)[i] + f(y0)[i]);
        }
    return ans;
}

// Матрица Якоби
vector<vector<TT>> Jacobi_matr(size_t const dim, const string &method, const vector<TT> &x0, const vector<TT> &y0) {
    vector<vector<TT>> Jacobi(dim, vector<TT>(dim));
    vector<TT> fx = calcMethod(dim, method, x0, y0);
    vector<TT> x0Vect(dim);
    vector<TT> xDiffVect(dim);

    for (size_t i = 0; i < dim; i++) {
        x0Vect = x0;
        x0Vect[i] += COMPARE_RATE;
        xDiffVect = calcMethod(dim, method, x0Vect, y0);
        for (size_t j = 0; j < dim; j++) {
            Jacobi[j][i] = (xDiffVect[j] - fx[j]) / COMPARE_RATE;
        }
    }
    return Jacobi;
}


// Метод Ньютона
vector<TT> Newton(const string &method, const size_t dim, const vector<TT> &y0) {
    vector<vector<TT>> Jacobi(dim, vector<TT>(dim));
    vector<TT> ans(dim);
    vector<TT> x(dim);
    vector<TT> xk(dim);

    do {
        x = xk;
        Jacobi = Jacobi_matr(dim, method, x, y0);

        if (dim == 2) {
            TT jacobian = Jacobi[0][0] * Jacobi[1][1] - Jacobi[1][0] * Jacobi[0][1];
            TT temp = Jacobi[0][0];
            Jacobi[0][0] = Jacobi[1][1] / jacobian;
            Jacobi[1][1] = temp / jacobian;
            Jacobi[0][1] /= -jacobian;
            Jacobi[1][0] /= -jacobian;
        } else Jacobi = inverseMatrix(Jacobi);

        xk = vectorMatrixMultiplication(Jacobi, calcMethod(dim, method, xk, y0));
        xk = vectorOperation(x, xk, '-');
    } while (norm1Vector(vectorOperation(x, xk, '-')) > COMPARE_RATE);
    ans = xk;
    return ans;
}
