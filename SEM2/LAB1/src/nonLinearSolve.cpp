#include <vector>
#include <string>
#include "shared.h"
#include "SOLE.h"

using namespace std;

vector<TT> f(const vector<TT>& x) {
    vector<TT> u(x.size());

    u[0] = x[1];
    u[1] = -20 / 0.3 * x[0];
    return u;
}

//template<typename T>
//using matrix = vector<vector<T>>;

/*auto derivative(function<TT(const vector<TT> &)> &fun, TT eps) {
    return [&fun, eps](vector<TT> &x) {
        vector<TT> xDiff(x);
        vectorDigit(eps, xDiff, '+');
        return ((fun(xDiff) - fun(x)) / eps);
    };
}

vector<TT> nonLinearSolve(function<TT(const vector<TT> &)> &fun,
                  const vector<TT> &initPoints, TT eps, size_t &iterCount) {
    vector<TT> currPoints = initPoints;
    auto diff = derivative(fun, COMPARE_RATE);
    TT diffVal;
    do {
        vector<TT> newPoints(currPoints);
        vectorDigit(fun(currPoints) / diff(currPoints), newPoints, '-');
        diffVal = normInfVector(vectorOperation(newPoints, currPoints, '-'));
        currPoints = newPoints;
        ++iterCount;
    } while (diffVal > eps);
    return currPoints;
}*/

vector<TT> calcMethod(const size_t dim, const string& method, const vector<TT>& x0, const vector<TT>& y0) {
    vector<TT> ans(dim);

    if (method == "Euler") {
        for (size_t i = 0; i < dim; i++) {
            ans[i] = x0[i] - y0[i] - step * f(x0)[i];
        }
    }

    if (method == "Symmetric")
        for (size_t i = 0; i < dim; i++)
        {
            ans[i] = x0[i] - y0[i] - step / 2 * (f(x0)[i] + f(y0)[i]);
        }
    return ans;
}

// Матрица Якоби
vector<vector<TT>> Jacobi_matr(size_t const dim, const string& method, const vector<TT>& x0, const vector<TT>& y0)
{
    TT eps1 = 1e-05;
    vector<vector<TT>> Jacobi(dim, vector<TT>(dim));
    vector<TT> fx = calcMethod(dim, method, x0, y0);
    vector<TT> temp1(dim);
    vector<TT> temp2(dim);

    for (size_t i = 0; i < dim; i++)
    {
        temp1 = x0;
        temp1[i] += eps1;
        temp2 = calcMethod(dim, method, temp1, y0);
        for (size_t j = 0; j < dim; j++)
        {
            Jacobi[j][i] = (temp2[j] - fx[j]) / eps1;
        }
    }
    return Jacobi;
}


// Метод Ньютона
vector<TT> Newton(const string& method, const size_t dim, const vector<TT>& y0)
{
    vector<vector<TT>> Jacobi(dim, vector<TT>(dim));
    vector<TT> ans(dim);
    vector<TT> x(dim);
    vector<TT> xk(dim);
    TT temp; TT jacobian;

    do
    {
        x = xk;
        Jacobi = Jacobi_matr(dim, method, x, y0);

        if (dim == 2)
        {
            jacobian = Jacobi[0][0] * Jacobi[1][1] - Jacobi[1][0] * Jacobi[0][1];
            temp = Jacobi[0][0];
            Jacobi[0][0] = Jacobi[1][1] / jacobian;
            Jacobi[1][1] = temp / jacobian;
            Jacobi[0][1] /= -jacobian;
            Jacobi[1][0] /= -jacobian;
        }
        else Jacobi = inverseMatrix(Jacobi);

        xk = vectorMatrixMultiplication(Jacobi, calcMethod(dim, method, xk, y0));
        xk = vectorOperation(x, xk, '-');
        x = vectorOperation(x, xk, '-');

        ans = xk;

    } while (norm1Vector(x) > COMPARE_RATE);
    return ans;
}
