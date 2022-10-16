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
    return currX;
}

vector<TT> setOmega(vector<vector<TT>> &matrix, vector<TT> &rightVect) {
    return relaxation(matrix, rightVect, 0.5);
}

vector<TT> zeidel(vector<vector<TT>> &matrix, vector<TT> &rightVect) {
    return relaxation(matrix, rightVect, 1);
}


void getMethod(vector<TT>(*func)(vector<vector<TT>> &, vector<TT> &), const string& method) {
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
