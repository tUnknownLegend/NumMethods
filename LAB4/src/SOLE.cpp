#include <iostream>
#include <string>
#include <cmath>
#include <<cassert>
#include "shared.h"

using std::vector;
using std::string;


void genHessenberg(vector<vector<TT>> &matrix, const unsigned int& k, const unsigned int& l){
    matrix = identityMatrix(matrix.size());
    TT alpha = matrix[k][k-1] / sqrt(pow(matrix[k][k-1], 2) + pow(matrix[l][k-1], 2));
    TT betha = matrix[l][k-1] / sqrt(pow(matrix[k][k-1], 2) + pow(matrix[l][k-1], 2));
    assert(k < l);
    assert(l < matrix.size());
    assert(std::abs(pow(alpha, 2) + pow(betha, 2)) <= COMPARE_RATE);
    matrix[k][k] = alpha;
    matrix[k][l] = betha;
    matrix[l][k] = -betha;
    matrix[l][l] = alpha;
}

vector<TT> Hessenberg(vector<vector<TT>> &matrix, vector<TT> &rightVect) {
    vector<vector<TT>> Tmatrix;
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            genHessenberg(Tmatrix, i, j);
            vector<vector<TT>> m = matrixOperations(Tmatrix, matrix, '*');
            matrix = matrixOperations(Tmatrix, transpoceMatrix(m), '*');
        }
    }
    return {};
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

void getHessenberg() {
    getMethod(&Hessenberg, "Hessenberg");
}