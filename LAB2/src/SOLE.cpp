#include <iostream>
#include "shared.h"

using std::vector;

const TT thau = 0.01;

vector<TT> iteration(vector<vector<TT>> &matrix, vector<TT> &rightVect) {
    matrixDigit(thau, matrix, '*');
    std::vector<std::vector<TT>> Cmatrix = matrixOperations(identityMatrix(rightVect.size()), matrix, '-');
    double norm = norm1Matrix(Cmatrix);
    vectorDigit(thau, rightVect, '*');

    vector<TT> prevX = rightVect;
    vector<TT> currX = vectorOperation(matrixVectorMultiplication(Cmatrix, prevX), rightVect, '+');
    if (norm <= 1) {
        while (norm1Vector(vectorOperation(currX, prevX, '-')) <= (1 - norm) / norm) {
            prevX = currX;
            currX = vectorOperation(matrixVectorMultiplication(Cmatrix, prevX), rightVect, '+');
        }
    } else {
        std::cout << "norm > 1\n";
    }
    return currX;
}

vector<TT> jacobi(vector<vector<TT>> &matrix, vector<TT> rightVect) {
    std::vector<std::vector<TT>> Cmatrix(matrix);
    vector<TT> Yvect(rightVect);

    for (int i = 0; i < rightVect.size(); ++i) {
        for (int j = 0; j < rightVect.size(); ++j) {
            Cmatrix[i][j] = -matrix[i][j] / matrix[i][i];
        }
        Yvect[i] = rightVect[i] / matrix[i][i];
    }
    matrixDigit(thau, Cmatrix, '*');
    double norm = norm1Matrix(Cmatrix);
    vector<TT> prevX(Yvect);
    vectorDigit(thau, prevX, '*');
    vector<TT> currX = vectorOperation(matrixVectorMultiplication(Cmatrix, prevX), rightVect, '+');
    if (norm <= 1) {
        while (norm1Vector(vectorOperation(currX, prevX, '-')) <= (1 - norm) / norm) {
            currX = vectorOperation(matrixVectorMultiplication(Cmatrix, prevX), rightVect, '+');
        }
    } else {
        std::cout << "norm > 1\n";
    }
    return currX;
}

void getMethod(vector<TT>(*func)(vector<vector<TT>> &, vector<TT> &)) {
    vector<vector<TT>> matrix;
    vector<TT> rightVect;
    inputMatrix(matrix);
    inputVector(rightVect);
    outputMatrix(matrix);
    rightVect = func(matrix, rightVect);
    outputVector(func(matrix, rightVect));
}

void getIterational() {
    getMethod(&iteration);
}

void getJacobi() {
    //getMethod(&jacobi);
}
