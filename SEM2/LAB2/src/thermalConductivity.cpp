#include <vector>
#include <iostream>

#include "shared.h"
#include "SOLE.h"

using namespace std;

TT u_0(const TT x) {
    return 1.0;
}

TT u_L(const TT x, const TT t) {
//    return 1.0;
    return sin(x) * exp(-t);
}

TT uInitial(const TT x, const TT L) {
    return 0.1 + x * (L - x);
}

// depends on coordinate
TT Kx(const TT x) {
    if (x <= x1) {
        return k1;
    } else if (x < x2) {
        return k1 * (x - x2) / (x1 - x2) + k2 * (x - x1) / (x2 - x1);
    } else {
        return k2;
    };
}

// depends on temperature
double Ku(const TT u, const TT kappa) {
    return kappa * pow(u, sigma);
}

double Pt(const TT t) {
    return 0.0;

    if (t >= 0 && t < t0) {
        return 2 * Q * t;
    } else {
        return 0;
    }
}

vector<TT>
tridiagonalMatrixAlgorithm(const int dim, const vector<TT> &a,
                           const vector<TT> &b, const vector<TT> &c_, const vector<TT> &d) {
    vector<TT> alpha(dim);
    vector<TT> beta(dim);

    alpha[0] = -c_[0] / b[0];
    beta[0] = d[0] / b[0];

    for (int i = 1; i < dim; i++) {
        alpha[i] = -c_[i] / (b[i] + a[i] * alpha[i - 1]);
        beta[i] = (-a[i] * beta[i - 1] + d[i]) / (a[i] * alpha[i - 1] + b[i]);
    }

    for (int i = dim - 2; i > -1; i--) {
        beta[i] = beta[i] + alpha[i] * beta[i + 1];
    }
    return beta;
}

typedef TT(*Fxt)(const TT x, const TT t);

typedef TT(*Fx)(const TT x);

auto doublePtMethod(Fx firstPt, Fx secondPt) {
    return [firstPt, secondPt](vector<vector<TT>> &Y, const vector<TT> &a,
                               const vector<TT> &diag1, vector<TT> &diag2, const vector<TT> &diag3,
                               const TT A0, const TT BN, const TT kappa, const TT secondKappa) {
        const TT L = N * h;
        diag2[0] += A0 * secondKappa;
        diag2[N - 2] += BN * kappa;

        vector<TT> right(N - 1);

        for (int j = 1; j < k; j++) {

            for (int i = 0; i < N - 1; i++)
                right[i] = -(c * ro * h / tao * Y[j - 1][i + 1] + (1 - sigma) * a[i] *
                                                                  (Y[j - 1][i + 2] -
                                                                   2 * Y[j - 1][i + 1] +
                                                                   Y[j - 1][i]) / h);

            TT mu1 = (c * ro * Y[j - 1][0] * h / (2 * tao) + sigma * firstPt(tao * j) +
                      (1 - sigma) * (firstPt(tao * (j - 1)) + (Y[j - 1][1] - Y[j - 1][0]) / h)) /
                     (c * ro * h / (2 * tao) + sigma * a[0] / h);
            TT mu2 = (c * ro * Y[j - 1][N] * h / (2 * tao) + sigma * secondPt(tao * j) +
                      (1 - sigma) * (secondPt(tao * (j - 1)) - (Y[j - 1][N] - Y[j - 1][N - 1]) / h)) /
                     (c * ro * h / (2 * tao) + sigma * a[N - 1] / h);

            right[0] -= A0 * mu1;
            right[N - 2] -= BN * mu2;

            vector<TT> temp = tridiagonalMatrixAlgorithm(N - 1, diag1, diag2, diag3, right);

            for (int i = 1; i < N; i++) {
                Y[j][i] = temp[i - 1];
            }

            Y[j][0] = kappa * Y[j][1] + mu1;
            Y[j][N] = secondKappa * Y[j][N - 1] + mu2;
        }
        return Y;
    };
}

auto singlePtMethod(Fxt u_L, Fx Pt) {
    return [u_L, Pt](vector<vector<TT>> &Y, const vector<TT> &a,
                     const vector<TT> &diag1, vector<TT> &diag2, const vector<TT> &diag3,
                     const TT A0, const TT BN, const TT kappa, const TT secondKappa) {
        const TT L = N * h;
        diag2[N - 2] += BN * kappa;

        vector<TT> right(N - 1);
        for (int j = 1; j < k; j++) {

            for (int i = 0; i < N - 1; i++) {
                right[i] = -(c * ro * h / tao * Y[j - 1][i + 1] +
                             (1 - sigma) * a[i] * (Y[j - 1][i + 2] - 2 * Y[j - 1][i + 1] + Y[j - 1][i]) / h);
            }

            const TT mu = (c * ro * Y[j - 1][0] * h / (2 * tao) + sigma * Pt(tao * j) +
                           (1 - sigma) * (Pt(tao * (j - 1)) + (Y[j - 1][1] - Y[j - 1][0]) / h)) /
                          (c * ro * h / (2 * tao) + sigma * a[0] / h);

            right[0] -= A0 * mu;
            right[N - 2] -= BN * u_L(L, j * tao);

            vector<TT> temp = tridiagonalMatrixAlgorithm(N - 1, diag1, diag2, diag3, right);

            for (int i = 1; i < N; i++)
                Y[j][i] = temp[i - 1];

            Y[j][0] = kappa * Y[j][1] + mu;
            Y[j][N] = u_L(L, j * tao);
        }
        return Y;
    };
}

auto constantTempMethod(Fx u_0, Fxt u_L) {
    return [u_0, u_L](vector<vector<TT>> &Y, const vector<TT> &a,
                      const vector<TT> &diag1, const vector<TT> &diag2, const vector<TT> &diag3,
                      const TT A0, const TT BN, const TT kappa, const TT secondKappa) {

        const TT L = N * h;

        vector<TT> right(N - 1);
        for (int j = 1; j < k; j++) {

            Y[j][0] = u_0(j * tao);

            for (int i = 0; i < N - 1; i++) {
                right[i] = -(c * ro * h / tao * Y[j - 1][i + 1] +
                             (1 - sigma) * a[i] * (Y[j - 1][i + 2] - 2 * Y[j - 1][i + 1] + Y[j - 1][i]) / h);
            }

            right[0] -= A0 * u_0(j * tao);
            right[N - 2] -= BN * u_L(L, j * tao);

            vector<TT> temp = tridiagonalMatrixAlgorithm(N - 1, diag1, diag2, diag3, right);

            for (int i = 1; i < N; i++)
                Y[j][i] = temp[i - 1];

            Y[j][N] = u_L(L, j * tao);
        }
        return Y;
    };
}

// Интегро-интерполяционный метод (для K = K(x))
vector<vector<TT>>
integroInterpolation(Fxt uInitial, Fx Kx, const auto &calcFunction) {
    vector<vector<TT>> Y(k, vector<TT>(N + 1));

    // начальные условия
    for (int i = 0; i <= N; i++) {
        Y[0][i] = uInitial(i * h, N * h);
    }

    vector<TT> a(N);
    // значения теплового потока в полуузлах
    for (int i = 0; i < N; i++) {
        a[i] = Kx((i + 0.5) * h);
    }

    const TT kappa = (sigma * a[0] / h) / (c * ro * h / (2 * tao) + sigma * a[0] / h);
    const TT secondKappa = (sigma * a[N - 1] / h) / (c * ro * h / (2 * tao) + sigma * a[N - 1] / h);

    vector<TT> diag1(N - 1);
    vector<TT> diag2(N - 1);
    vector<TT> diag3(N - 1);
    for (int i = 0; i < N - 1; i++) {
        diag1[i] = sigma / h * a[i];
        diag3[i] = sigma / h * a[i + 1];
        diag2[i] = -(diag1[i] + diag3[i] + c * ro * h / tao);
    }

    TT A0 = diag1[0];
    TT BN = diag3[N - 2];

    diag1[0] = 0.0;
    diag3[N - 2] = 0.0;

    return calcFunction(Y, a, diag1, diag2, diag3, A0, BN, kappa, secondKappa);
}

TT calcDiff(const vector<vector<TT>> &answer) {
    TT error = 0.0;
    for (int j = 0; j < k; ++j) {
        for (int i = 0; i < N + 1; ++i) {
            error = max(abs(answer[j][i] - (cos(i * h) + sin(i * h)) * exp(-j * tao)), error);
        }
    }
    return error;
}

void IntegroInterpolation() {
    const auto doublePtFunction =
            doublePtMethod(Pt, Pt);

    const auto singlePtFunction =
            singlePtMethod(u_L, Pt);

    const auto constantTempCalcFunction =
            constantTempMethod(u_0, u_L);

    const auto ansZero = integroInterpolation(uInitial, Kx, constantTempCalcFunction);
    const auto ansSingle = integroInterpolation(uInitial, Kx, singlePtFunction);
    const auto ansDouble = integroInterpolation(uInitial, Kx, doublePtFunction);

    std::cout << calcDiff(ansZero) << "\n";
    std::cout << calcDiff(ansSingle) << "\n";
    std::cout << calcDiff(ansDouble) << "\n";

    outputMatrix(ansZero, "../data/zeroPt.txt");
    outputMatrix(ansSingle, "../data/singlePt.txt");
    outputMatrix(ansDouble, "../data/doublePt.txt");
}