#include <iostream>
#include <string>
#include <cmath>
#include <cassert>
#include "shared.h"

using std::vector;
using std::string;

void getNextMatrix(vector<vector<TT>>& matrix, const bool& last = false) {
    vector<vector<TT>> R(matrix);
    vector<vector<TT>> T = identityMatrix(matrix.size());
    vector<TT> rightVect(matrix.size(), 0);
    TT c;
    TT s;
    vector<TT> resultVect(rightVect.size(), 0.0);
    vector<TT> add(rightVect);

    for (int i = 0; i < rightVect.size(); ++i) {  // - 1
        for (int j = i + 1; j < rightVect.size(); ++j) {
            if (std::abs(R[i][j]) > COMPARE_RATE) {
                TT tempres = sqrt(pow(R[i][i], 2) + pow(R[j][i], 2));
                c = R[i][i] / tempres;
                s = R[j][i] / tempres;
                for (int k = 0; k < rightVect.size(); ++k) {
                    TT tr = R[i][k];
                    R[i][k] = c * R[i][k] + s * R[j][k];
                    R[j][k] = -s * tr + c * R[j][k];

                    TT tq = T[i][k];
                    T[i][k] = c * T[i][k] + s * T[j][k];
                    T[j][k] = -s * tq + c * T[j][k];
                }
                TT tv = add[i];
                add[i] = c * add[i] + s * add[j];
                add[j] = -s * tv + c * add[j];
            }
        }
    }

    for (int i = rightVect.size() - 1; i >= 0; --i) {
        TT sumJ = 0.0;
        for (int j = i + 1; j < rightVect.size(); ++j) {
            sumJ += R[i][j] * resultVect[j];
        }
        resultVect[i] = (add[i] - sumJ) / R[i][i];
    }

    for (int i = 0; i < matrix.size(); ++i) {
        if (std::abs(R[i][i]) < COMPARE_RATE) {
            std::cerr << "Matrix is singular";
            return;
        }
    }
    vector<vector<TT>> Q = transpoceMatrix(T);

//    std::cout << "Q:" << std::endl;
//    outputOnTheScreenMatrix(Q);
//    std::cout << "R:" << std::endl;
//    outputOnTheScreenMatrix(R);

    //std::cout << "\nAmount of mult: " << coounter << "\n";
    //std::cout << "\nAmount of mltB: " << mltB << "\n";

    matrix = (last ? R : matrixOperations(R, Q, '*'));
}

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

vector<TT> qrEigenvalues(vector<vector<TT>>& matrix, vector<TT> & vec) {
    for (int i = 0; i < matrix.size(); ++i) {
        getNextMatrix(matrix);
    }
    // outputOnTheScreenMatrix(matrix);

    vector<TT> resVec(matrix.size(), 0);

    for (int i = 0; i < matrix.size(); ++i){
            resVec [i] = matrix[i][i];
    }

    return resVec;
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
    // getMethod(&Hessenberg, "Hessenberg");
}

void getEigenvalues() {
    getMethod(&qrEigenvalues, "qrEigenvalues");
}
