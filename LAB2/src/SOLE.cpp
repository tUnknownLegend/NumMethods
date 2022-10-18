#include <iostream>
#include <string>
#include <cmath>
#include "shared.h"

using std::vector;
using std::string;

const TT thau = 10e-3;
const TT eps = 10e-4;
const TT eps0 = 10e-4;
//const vector<TT> ans{5, -7, 12, 4};
const vector<TT> ans{10, -10, 12, 4};

vector<TT> iteration(vector<vector<TT>> &matrix, vector<TT> &rightVect) {
    matrixDigit(thau, matrix, '*');
    std::vector<std::vector<TT>> Cmatrix = matrixOperations(identityMatrix(rightVect.size()), matrix, '-');
    vectorDigit(thau, rightVect, '*');

    //double norm = normInfMatrix(Cmatrix);
    double norm = norm1Matrix(Cmatrix);
    std::cout << "norm of matrix C: " << norm << "\n";
    vector<TT> prevX = rightVect;
    vector<TT> currX = vectorOperation(matrixVectorMultiplication(Cmatrix, prevX), rightVect, '+');
    vector<TT> initX = currX;
    if (norm < 1) {
        unsigned int amtIt = 0;
        // 1: norm1Vector(vectorOperation(currX, prevX, '-')) > norm1Vector(currX)
        // 2: norm1Vector(vectorOperation(currX, prevX, '-')) > norm1Vector(currX) * eps + eps0
        // 3: norm1Vector(vectorOperation(currX, prevX, '-')) > (1 - norm) / norm * eps
        while (norm1Vector(vectorOperation(currX, prevX, '-')) > ((1 - norm) / norm) * eps) {
            ++amtIt;
            prevX = currX;
            currX = vectorOperation(matrixVectorMultiplication(Cmatrix, prevX), rightVect, '+');
        }
        std::cout << "iterations: " << amtIt << "\n";
    } else {
        std::cerr << "norm > 1\n";
    }

    std::cout << "kEst: " << log((1 - norm) * eps / norm1Vector(vectorOperation(currX, initX, '-'))) / log(norm)
              << "\n";
    std::cout << "norm of error: " << norm1Vector(vectorOperation(ans, currX, '-')) << "\n";

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
    std::cout << "norm of matrix C: " << norm << "\n";
    vector<TT> prevX = rightVect;
    vector<TT> currX = vectorOperation(matrixVectorMultiplication(Cmatrix, prevX), Yvect, '+');
    vector<TT> initX = currX;
    if (norm < 1) {
        unsigned int amtIt = 0;
        // 1: norm1Vector(vectorOperation(currX, prevX, '-')) > norm1Vector(currX)
        // 2: norm1Vector(vectorOperation(currX, prevX, '-')) > norm1Vector(currX) * eps + eps0
        // 3: norm1Vector(vectorOperation(currX, prevX, '-')) > (1 - norm) / norm * eps
        while (norm1Vector(vectorOperation(currX, prevX, '-')) > (1 - norm) / norm * eps) {
            ++amtIt;
            prevX = currX;
            currX = vectorOperation(matrixVectorMultiplication(Cmatrix, prevX), Yvect, '+');
        }
        std::cout << "iterations: " << amtIt << "\n";
    } else {
        std::cerr << "norm > 1\n";
    }

    std::cout << "kEst: " << log((1 - norm) * eps / norm1Vector(vectorOperation(currX, initX, '-'))) / log(norm)
              << "\n";
    std::cout << "norm of error: " << norm1Vector(vectorOperation(ans, currX, '-')) << "\n";
    return currX;
}

vector<TT> getInitX(const vector<vector<TT>> &matrix, const vector<TT> &rightVect, const TT &omega) {
    vector<TT> currX(rightVect.size(), 0);
    vector<TT> prevX(rightVect.size(), 1);
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
    return currX;
}

vector<TT> relaxation(vector<vector<TT>> &matrix, vector<TT> &rightVect, const TT &omega) {
    vector<TT> currX(rightVect.size(), 0);
    vector<TT> prevX(rightVect.size(), 1);

    //  geting C matrix
    vector<vector<TT>> Lmatrix(matrix);
    vector<vector<TT>> Dmatrix(matrix);
    vector<vector<TT>> Umatrix(matrix);
    LDU(matrix, Lmatrix, Dmatrix, Umatrix);
    vector<vector<TT>> inverseD = inverseMatrix(Dmatrix);
    vector<vector<TT>> E = identityMatrix(matrix.size());

    vector<TT> initX = getInitX(matrix, rightVect, omega);

    Dmatrix = matrixOperations(inverseD, Lmatrix, '*');
    matrixDigit(omega, Dmatrix, '*');
    vector<vector<TT>> Cl = matrixOperations(E, Dmatrix, '+');

    Dmatrix = matrixOperations(inverseD, Umatrix, '*');
    matrixDigit(omega, Dmatrix, '*');
    if (omega == 1) {
        matrixDigit((1 - omega), E, '*');
    }
    vector<vector<TT>> Cu = matrixOperations(E, Dmatrix, '-');

    TT norm = norm1Matrix(matrixOperations(inverseMatrix(Cl), Cu, '*'));
    std::cout << "norm of matrix C: " << norm << "\n";
/*    std::cout << "Cmatrix: ";
    outputOnTheScreenMatrix(Cmatrix);
    std::cout << "Cu: ";
    outputOnTheScreenMatrix(Cu);
    std::cout << "\n";*/

    unsigned int amtIt = 0;

    // 1: norm1Vector(vectorOperation(currX, prevX, '-')) > norm1Vector(currX)
    // 2: norm1Vector(vectorOperation(currX, prevX, '-')) > norm1Vector(currX) * eps + eps0
    // 3: norm1Vector(vectorOperation(currX, prevX, '-')) > (1 - norm) / norm1Matrix(Cu) * eps
    while (norm1Vector(vectorOperation(currX, prevX, '-')) > (1 - norm) / norm1Matrix(Cu) * eps) {
        prevX = currX;
        ++amtIt;
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
    std::cout << "iterations: " << amtIt << "\n";
    std::cout << "kEst: " << log((1 - norm) * eps / norm1Vector(vectorOperation(currX, initX, '-'))) / log(norm)
              << "\n";
    std::cout << "norm of error: " << norm1Vector(vectorOperation(ans, currX, '-')) << "\n";
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

    unsigned int amtIt = 0;
    while (norm1Vector(vectorOperation(currX, prevX, '-')) > norm1Vector(currX) * eps + eps0) {
        prevX = currX;
        ++amtIt;
        currX[0] = (1 - omega) * prevX[0] + omega * (-rightDiag[0] / midDiag[0] * currX[1] + vect[0] / midDiag[0]);
        for (int i = 1; i < vect.size() - 1; ++i) {
            currX[i] = omega * (-leftDiag[i] / midDiag[i] * currX[i - 1] -
                                rightDiag[i] / midDiag[i] * currX[i + 1] + vect[i] / midDiag[i]) +
                       (1 - omega) * prevX[i];
        }
        int last = vect.size() - 1;
        currX[last] = omega * (-leftDiag[last] / midDiag[last] * currX[last - 1] + vect[last] / midDiag[last]) +
                      (1 - omega) * prevX[last];
    }
    std::cout << "iterations: " << amtIt << "\n";
    //  getting C matrix

    vector<TT> Cmatrix(leftDiag.size(), 0);
    vector<TT> Cu(leftDiag.size(), 0);
    for (int i = 0; i < midDiag.size(); ++i) {
        Cmatrix[i] = -omega / midDiag[i] * leftDiag[i];
        Cu[i] = -omega / midDiag[i] * rightDiag[i];
    }


    std::cout << "norm of matrix Cmatrix: " << norm1Vector(Cmatrix) << "\n";
    std::cout << "norm of matrix Cu: " << norm1Vector(Cu) << "\n";
//    std::cout << "Cmatrix: ";
//    outputOnTheScreenVector(Cmatrix);
//    std::cout << "Cu: ";
//    outputOnTheScreenVector(Cu);
//    std::cout << "\n";

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

    std::cout << "three diag:\n";
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
    std::cout << "----------------------\n";
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
