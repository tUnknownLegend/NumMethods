#include <iostream>
#include <string>
#include "shared.h"

using std::vector;
using std::string;

const TT thau = 10e-3;
const TT eps = 10e-5;
const TT eps0 = 10e-5;

vector<TT> iteration(vector<vector<TT>> &matrix, vector<TT> &rightVect) {
    matrixDigit(thau, matrix, '*');
    std::vector<std::vector<TT>> Cmatrix = matrixOperations(identityMatrix(rightVect.size()), matrix, '-');
    vectorDigit(thau, rightVect, '*');

    //double norm = normInfMatrix(Cmatrix);
    double norm = norm1Matrix(Cmatrix);
    vector<TT> prevX = rightVect;
    vector<TT> currX = vectorOperation(matrixVectorMultiplication(Cmatrix, prevX), rightVect, '+');
    if (norm < 1) {
        while (norm1Vector(vectorOperation(currX, prevX, '-')) > ((1 - norm) / norm) * eps) {
            prevX = currX;
            currX = vectorOperation(matrixVectorMultiplication(Cmatrix, prevX), rightVect, '+');
        }
    } else {
        std::cerr << "norm > 1\n";
    }
    return currX;
}

vector<TT> jacobi(vector<vector<TT>> &matrix, vector<TT> &rightVect) {
    std::vector<std::vector<TT>> Cmatrix(matrix);
    vector<TT> Yvect(rightVect);

    for (int i = 0; i < rightVect.size(); ++i) {
        for (int j = 0; j < rightVect.size(); ++j) {
            if (i == j) {
                Cmatrix[i][j] = 0;
            } else {
                Cmatrix[i][j] = -matrix[i][j] / matrix[i][i];
            }
        }
        Yvect[i] = rightVect[i] / matrix[i][i];
    }

    //double norm = normInfMatrix(Cmatrix);
    double norm = norm1Matrix(Cmatrix);
    vector<TT> prevX = rightVect;
    vector<TT> currX = vectorOperation(matrixVectorMultiplication(Cmatrix, prevX), Yvect, '+');
    if (norm < 1) {
        while (norm1Vector(vectorOperation(currX, prevX, '-')) > (1 - norm) / norm * eps) {
            prevX = currX;
            currX = vectorOperation(matrixVectorMultiplication(Cmatrix, prevX), Yvect, '+');
        }
    } else {
        std::cerr << "norm > 1\n";
    }
    return currX;
}

vector<TT> relaxation(vector<vector<TT>> &matrix, vector<TT> &rightVect, const TT &omega) {
    vector<TT> currX(rightVect.size(), 0);
    vector<TT> prevX(rightVect.size(), 1);

    while (norm1Vector(vectorOperation(currX, prevX, '-')) > norm1Vector(currX) * eps + eps0) {
        prevX = currX;
        for (int i = 0; i < rightVect.size(); ++i) {
            TT prevSum = 0.0;
            TT currSum = 0.0;

            for (int j = i + 1; j < rightVect.size(); ++j) {
                prevSum += matrix[i][j] * prevX[j];
            }
            for (int j = 0; j < i; ++j) {
                currSum += matrix[i][j] * currX[j];
            }
            currX[i] =
                    -omega * (currSum + prevSum - rightVect[i]) / matrix[i][i] + (1 - omega) * prevX[i];
        }
    }
    //  geting C matrix
    vector<vector<TT>> Lmatrix(matrix);
    vector<vector<TT>> Dmatrix(matrix);
    vector<vector<TT>> Umatrix(matrix);
    LDU(matrix, Lmatrix, Dmatrix, Umatrix);
    vector<vector<TT>> inverseD = inverseMatrix(Dmatrix);
    std::move(Dmatrix);

    vector<vector<TT>> Cmatrix(matrix);
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            Cmatrix[i][j] = -omega * inverseD[i][i] * Lmatrix[i][j];
        }
    }

    vector<vector<TT>> Cu(matrix);
    vector<vector<TT>> E = identityMatrix(matrix.size());
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = matrix.size() - 1; j > i; --j) {
            if (i == j) {
                Cu[i][j] = (1 - omega) - omega * inverseD[i][i] * Umatrix[i][j];
            } else {
                Cu[i][j] = -omega * inverseD[i][i] * Umatrix[i][j];
            }
        }
    }

    //outputMatrix(Cmatrix);
    //outputMatrix(Cu);
    return currX;
}

vector<TT> setOmega(vector<vector<TT>> &matrix, vector<TT> &rightVect) {
    return relaxation(matrix, rightVect, 0.5);
}

vector<TT> zeidel(vector<vector<TT>> &matrix, vector<TT> &rightVect) {
    return relaxation(matrix, rightVect, 1);
}


vector<TT>
threeDiag(vector<TT> &leftDiag, vector<TT> &midDiag, vector<TT> &rightDiag,
          vector<TT> &vect, const TT &omega) {
    vector<TT> currX(vect.size(), 0);
    vector<TT> prevX(vect.size(), 1);

    while (norm1Vector(vectorOperation(currX, prevX, '-')) > norm1Vector(currX) * eps + eps0) {
        prevX = currX;

        currX[0] = (1 - omega) * prevX[0] + omega * (-rightDiag[0] / midDiag[0] * currX[1] + vect[0] / midDiag[0]);
        for (int i = 1; i < vect.size() - 1; ++i) {
            currX[i] = omega * (-leftDiag[i] / midDiag[i] * currX[i - 1] -
                                rightDiag[i] / midDiag[i] * currX[i + 1] + vect[i] / midDiag[i]) +
                       (1 - omega) * prevX[i];
        }
        int last = vect.size() - 1;
        currX[last] = omega * (-leftDiag[last] / midDiag[last] * currX[last - 1] + vect[last] / midDiag[last]) + (1 - omega) * prevX[last];
    }
    //  getting C matrix

    vector<TT> Cmatrix(leftDiag.size(), 0);
    vector<TT> Cu(leftDiag.size(), 0);
    for (int i = 0; i < midDiag.size(); ++i) {
        Cmatrix[i] = -omega / midDiag[i] * leftDiag[i];
        Cu[i] = -omega / midDiag[i] * rightDiag[i];
    }

    return currX;
}

void getThreeDiag() {
    const int size = 4;
    vector<TT> leftDiag(size, 1.0);
    vector<TT> midDiag(size, 4.0);
    vector<TT> rightDiag(size, 1.0);
    vector<TT> vect(size, 6.0);

    for (int i = 0; i < vect.size(); ++i) {
        vect[i] = 10 - 2 * ((i + 1) % 2);
    }
    vect[0] = 6;
    vect[vect.size() - 1] = 9 - 3 * (vect.size() % 2);

    leftDiag[0] = 0;
    rightDiag[rightDiag.size() - 1] = 0;

    std::cout << "three diag: ";
    outputOnTheScreenVector(threeDiag(leftDiag, midDiag, rightDiag, vect, 1.0));
}

void getMethod(vector<TT>(*func)(vector<vector<TT>> &, vector<TT> &), const string &method) {
    vector<vector<TT>> matrix;
    vector<TT> rightVect;
    inputMatrix(matrix);
    inputVector(rightVect);
    outputMatrix(matrix);
    vector<TT> res = func(matrix, rightVect);
    std::cout << method << ": ";
    outputOnTheScreenVector(res);
    outputVector(res);
}

void getIterational() {
    getMethod(&iteration, "iteration");
}

void getJacobi() {
    getMethod(&jacobi, "jacobi");
}

void getRelaxation() {
    getMethod(&setOmega, "relaxation");
}

void getZeidel() {
    getMethod(&zeidel, "zeidel");
}
