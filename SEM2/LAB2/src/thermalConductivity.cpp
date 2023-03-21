#include <cmath>
#include <vector>

#include "shared.h"
#include "SOLE.h"

using namespace std;

TT u_x0(const TT x) {
    return 1.0;
}

TT u_Lt(const TT x, const TT t) {
    return 1.0;
}

TT ux0_0(const TT x, const TT u0) {
    return 1 - x * (1 - x);
    //return 0.0;
}

TT Kx(const TT x) {
    const TT k1 = 0.1;
    const TT k2 = 1;
    const TT x1 = 0.0;
    const TT x2 = 1.0;
    const TT L = 1;

    //if (x >= 0 && x <= L)
    //{
    //	if (x <= x1)
    //	{
    //		return k1;
    //	}
    //	else if (x < x2)
    //	{
    //		return k1 * (x - x2) / (x1 - x2) + k2 * (x - x1) / (x2 - x1);
    //	}
    //	else return k2;
    //}
    //else return 0;
    return 4;
}

double Pt(const TT t) {
    const TT t0 = 0.5;
    const TT Q = 10;
    //if (t >= 0 && t < t0)
    //{
    //	//return 0;
    //	return 2 * Q * t;
    //}
    //else return 0;
    return 0.0;
}

vector<TT>
tridiagonalMatrixAlgorithm(int dim, vector<TT> a, vector<TT> b, vector<TT> c, vector<TT> d) {
    vector<TT> alpha(dim);
    vector<TT> beta(dim);

    alpha[0] = -c[0] / b[0];
    beta[0] = d[0] / b[0];

    for (int i = 1; i < dim; i++) {
        alpha[i] = -c[i] / (b[i] + a[i] * alpha[i - 1]);
        beta[i] = (-a[i] * beta[i - 1] + d[i]) / (a[i] * alpha[i - 1] + b[i]);
    }

    for (int i = dim - 2; i > -1; i--) {
        beta[i] = beta[i] + alpha[i] * beta[i + 1];
    }
    return beta;
}

typedef TT(*Fxt)(const TT x, const TT t);

typedef TT(*Fx)(const TT x);

auto defaultMethod(const int k, const TT h, const TT tao, const TT c, const TT p, Fxt uLt, Fx Pt) {
    return [k, h, tao, uLt, c, p, Pt](int N,
                                      vector<vector<TT>> &Y, const vector<TT> &a,
                                      const vector<TT> &diag1, const vector<TT> &diag2, const vector<TT> &diag3,
                                      const TT A0, const TT BN, const TT kappa, vector<TT> &right) {
        const TT L = N * h;

        for (int j = 1; j < k; j++) {

            for (int i = 0; i < N - 1; i++) {
                right[i] = -(c * p * h / tao * Y[j - 1][i + 1] +
                             (1 - sigma) * a[i] * (Y[j - 1][i + 2] - 2 * Y[j - 1][i + 1] + Y[j - 1][i]) / h);
            }

            TT mu = (c * p * Y[j - 1][0] * h / (2 * tao) + sigma * Pt(tao * j) +
                     (1 - sigma) * (Pt(tao * (j - 1)) + (Y[j - 1][1] - Y[j - 1][0]) / h)) /
                    (c * p * h / (2 * tao) + sigma * a[0] / h);

            right[0] -= A0 * mu;
            right[N - 2] -= BN * uLt(L, j * tao);

            vector<TT> temp = tridiagonalMatrixAlgorithm(N - 1, diag1, diag2, diag3, right);

            for (int i = 1; i < N; i++)
                Y[j][i] = temp[i - 1];

            Y[j][0] = kappa * Y[j][1] + mu;
            Y[j][N] = uLt(L, j * tao);
        }
        return Y;
    };
}

auto constantTempMethod(int k, const TT h, const TT tao, const TT c, const TT p, Fx u_0t) {
    return [k, h, tao, c, p, u_0t](int N,
                                   vector<vector<TT>> &Y, const vector<TT> &a,
                                   const vector<TT> &diag1, const vector<TT> &diag2, const vector<TT> &diag3,
                                   const TT A0, const TT BN, const TT kappa, vector<TT> &right) {

        const TT L = N * h;
        for (int j = 1; j < k; j++) {

            Y[j][0] = u_0t(j * tao);

            for (int i = 0; i < N - 1; i++) {
                right[i] = -(c * p * h / tao * Y[j - 1][i + 1] +
                             (1 - sigma) * a[i] * (Y[j - 1][i + 2] - 2 * Y[j - 1][i + 1] + Y[j - 1][i]) / h);
            }

            right[0] -= A0 * u_0t(j * tao);
            right[N - 2] -= BN * u_Lt(L, j * tao);

            vector<TT> temp = tridiagonalMatrixAlgorithm(N - 1, diag1, diag2, diag3, right);

            for (int i = 1; i < N; i++)
                Y[j][i] = temp[i - 1];

            Y[j][N] = u_Lt(L, j * tao);
        }
        return Y;
    };
}

// Интегро-интерполяционный метод (для K = K(x))
vector<vector<TT>>
integroInterpolation(const int n, const int k, const TT h, const TT tao, const TT c, const TT p, Fxt ux0, Fx Kx,
                     const auto &calcFunction) {
    const int N = n - 1;
    vector<vector<TT>> Y(k, vector<TT>(n));

    // начальные условия
    for (int i = 0; i < n; i++) {
        Y[0][i] = ux0(i * h, (n - 1) * h);
    }

    vector<TT> a(N);
    // значения теплового потока в полуузлах
    for (int i = 0; i < N; i++) {
        a[i] = Kx((i + 0.5) * h);
    }

    const TT kappa = (sigma * a[0] / h) / (c * p * h / (2 * tao) + sigma * a[0] / h);

    vector<TT> diag1(N - 1);
    vector<TT> diag2(N - 1);
    vector<TT> diag3(N - 1);
    for (int i = 0; i < N - 1; i++) {
        diag1[i] = sigma / h * a[i];
        diag3[i] = sigma / h * a[i + 1];
        diag2[i] = -(diag1[i] + diag3[i] + c * p * h / tao);
    }

    TT A0 = diag1[0];
    TT BN = diag3[N - 2];

    diag1[0] = 0.0;
    diag3[N - 2] = 0.0;
    diag2[N - 2] += BN * kappa;

    vector<TT> right(N - 1);

    return calcFunction(N, Y, a, diag1, diag2, diag3, A0, BN, kappa, right);
}

void IntegroInterpolation(const bool isDefault) {
    TT ro = 4;
    TT c = 0.5;
    TT t0 = 0.5;
    TT l = 1.0;

    TT h = 0.08;
    TT tao = 0.002;
    int n = round(l / h + 1);
    int k = round(t0 / tao + 1);

    const auto defaultCalcFunction =
            defaultMethod(k, h, tao, c, ro, u_Lt, Pt);

    const auto constantTempCalcFunction =
            constantTempMethod(k, h, tao, c, ro, u_x0);

    const auto answer = (isDefault ?
                         integroInterpolation(n, k, h, tao, c, ro, ux0_0, Kx, defaultCalcFunction) :
                         integroInterpolation(n, k, h, tao, c, ro, ux0_0, Kx, constantTempCalcFunction));

    (isDefault ? outputMatrix(answer, "../defaultMethod.txt") :
     outputMatrix(answer, "../constTempMethod.txt"));
}