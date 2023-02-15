#include <cmath>
#include "shared.h"
#include "nonLinearSolve.h"
#include "SOLE.h"

using namespace std;

vector<vector<TT>> explicitEuler(const vector<TT> &cond, int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    for (int i = 0; i < n - 1; ++i) {
        y[i + 1] = vectorOperation(y[i],
                                     vectorRDigit(step, f(y[i]), '*'), '+');
    }
    return y;
}

vector<vector<TT>> implicitEuler(const vector<TT> &cond, const int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    for (int i = 0; i < n - 1; i++) {
        y[i + 1] = Newton("Euler", cond.size(), y[i]);
    }
    return y;
}

vector<vector<TT>> symmetric(const vector<TT> &cond, const int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    for (int i = 0; i < n - 1; i++) {
        y[i + 1] = Newton("Symmetric", cond.size(), y[i]);
    }
    return y;
}

void k_iSum(vector<TT>& temp, const vector<TT>& k1, const vector<TT>& k2, const vector<TT>& k3,
            const vector<TT>& k4, const vector<TT>& y_i) {
    temp = vectorOperation(vectorRDigit(2.0, k2, '*'), vectorRDigit(2.0, k3, '*'), '+');
    temp = vectorOperation(k1, temp, '+');
    temp = vectorOperation(k4, temp, '+');
    temp = vectorRDigit(step / 6.0, temp, '*');
    temp = vectorOperation(y_i, temp, '+');
}

vector<vector<TT>> rungeKutta(const vector<TT> &cond, const int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    vector<TT> k1(cond.size());
    vector<TT> k2(cond.size());
    vector<TT> k3(cond.size());
    vector<TT> k4(cond.size());
    vector<TT> temp1(cond.size());
    vector<TT> temp2(cond.size());
    const TT tau = 1.0;
    for (int i = 0; i < n - 1; i++) {
        // calc k_i
        k1 = f(y[i]);
        temp1 = vectorOperation(y[i], vectorRDigit(0.5 * step, k1, '*'), '+');
        k2 = f(temp1);
        temp1 = vectorOperation(y[i], vectorRDigit(0.5 * step, k2, '*'), '+');
        k3 = f(temp1);
        temp1 = vectorOperation(y[i], vectorRDigit(1.0 * step, k3, '*'), '+');
        k4 = f(temp1);

        // sum k_i
        k_iSum(temp1, k1, k2, k3, k4, y[i]);

        temp2 = y[i];
        for (int j = 0; j < 2; j++) {
            k1 = f(y[i]);
            temp2 = vectorOperation(temp2,
                                    vectorRDigit(tau / 4 * step, k1, '*'), '+');
            k2 = f(temp2);
            temp2 = vectorOperation(temp2,
                                    vectorRDigit(tau / 4 * step, k2, '*'), '+');
            k3 = f(temp2);
            temp2 = vectorOperation(temp2,
                                    vectorRDigit(tau / 2 * step, k3, '*'), '+');
            k4 = f(temp2);

            // k_i sum
            k_iSum(temp2, k1, k2, k3, k4, y[i]);
        }
        if (norm1Vector(vectorRDigit(1 / (pow(2, 4) - 1),
                                     (vectorOperation(temp2, temp1, '-')), '*'))
            <= COMPARE_RATE) {
            y[i + 1] = temp2;
        } else {
            y[i + 1] = temp1;
        }
    }
    return y;
}

void ImplicitEuler() {
    outputOnTheScreenMatrix(explicitEuler(initPoints, sizeN));
}

void ExplicitEuler() {
    outputOnTheScreenMatrix(implicitEuler(initPoints, sizeN));
}

void Symmetric() {
    outputOnTheScreenMatrix(symmetric(initPoints, sizeN));
}

void RungeKutta() {
    outputOnTheScreenMatrix(rungeKutta(initPoints, sizeN));
}

