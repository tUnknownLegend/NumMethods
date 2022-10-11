#include <iostream>
#include "shared.h"

using std::vector;

const TT thau = 0.01;

vector<TT> iteration(const vector<vector<TT>> &matrix, const vector<TT> &rightVect) {
    std::vector<std::vector<TT>> Ematrix = identityMatrix(rightVect.size());
    std::vector<std::vector<TT>> Cmatrix = matrixOperations(Ematrix, matrix, '-');
    matrixDigit(thau, Cmatrix, '*');
    double norm = norm1Matrix(Cmatrix);
    vector<TT> prevX(rightVect);
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

void getIterational() {
    vector<vector<TT>> matrix;
    vector<TT> rightVect;
    inputMatrix(matrix);
    inputVector(rightVect);
    outputMatrix(matrix);
    //outputVector(rightVect);

    outputVector(iteration(matrix, rightVect));
}