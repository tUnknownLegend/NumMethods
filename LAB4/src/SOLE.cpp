﻿#include <iostream>
#include <string>
#include <cmath>
#include <cassert>
#include "shared.h"

using std::vector;
using std::string;

const TT theta = 5.0;

std::pair<vector<vector<TT>>, vector<vector<TT>>>QR(const vector<vector<TT>>& matrix) {
    vector<vector<TT>> R(matrix);
    vector<vector<TT>> T = identityMatrix(matrix.size());
    vector<TT> rightVect(matrix.size(), 0);
    TT c;
    TT s;
    vector<TT> resultVect(rightVect.size(), 0.0);
    vector<TT> add(rightVect);

    for (int i = 0; i < matrix.size(); ++i) {  // - 1
        for (int j = i + 1; j < matrix.size(); ++j) {
            if (std::abs(R[i][j]) > COMPARE_RATE) {
                TT tempres = sqrt(pow(R[i][i], 2) + pow(R[j][i], 2));
                c = R[i][i] / tempres;
                s = R[j][i] / tempres;
                for (int k = 0; k < matrix.size(); ++k) {
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

    for (int i = matrix.size() - 1; i >= 0; --i) {
        TT sumJ = 0.0;
        for (int j = i + 1; j < matrix.size(); ++j) {
            sumJ += R[i][j] * resultVect[j];
        }
        resultVect[i] = (add[i] - sumJ) / R[i][i];
    }

    for (int i = 0; i < matrix.size(); ++i) {
        if (std::abs(R[i][i]) < COMPARE_RATE) {
            std::cerr << "Matrix is singular";
            return {};
        }
    }
    vector<vector<TT>> Q = transpoceMatrix(T);

    return {Q, R};
}

void getNextMatrix(vector<vector<TT>>& matrix, const bool& last = false) {
    vector<vector<TT>>tempMatrix(matrix);
    vector<vector<TT>> thetaMatrix = identityMatrix(
            matrix.size(), theta);
    matrix = matrixOperations(tempMatrix, thetaMatrix, '-');

    std::pair<vector<vector<TT>>, vector<vector<TT>>> qr = QR(matrix);

    tempMatrix = (last ? qr.second : matrixOperations(qr.second, qr.first, '*'));
    matrix = (last ? tempMatrix : matrixOperations(tempMatrix, thetaMatrix, '+'));
}

void genHessenberg(vector<vector<TT>> &matrix, const unsigned int& k, const unsigned int& l){
    TT temp = sqrt(pow(matrix[k][k-1], 2) + pow(matrix[l][k-1], 2));
    TT alpha = matrix[k][k-1] / temp;
    TT betha = matrix[l][k-1] / temp;
//    assert(std::abs(pow(alpha, 2) + pow(betha, 2) - 1) <= DIVISTION_ERROR);

    for (int i = k - 1; i < matrix.size(); ++i) {
        TT Mki = matrix[k][i];
        matrix[k][i] = alpha * matrix[k][i] + betha * matrix[l][i];
        matrix[l][i] = alpha * matrix[l][i] - betha * Mki;
    }
    for (int i = k - 1; i < matrix.size(); ++i) {
        TT Mik = matrix[i][k];
        matrix[i][k] = alpha * matrix[i][k] + betha * matrix[i][l];
        matrix[i][l] = alpha * matrix[i][l] - betha * Mik;
    }
}

void Hessenberg(vector<vector<TT>> &matrix) {
    for (int i = 1; i < matrix.size(); ++i) {
        for (int j = i + 1; j < matrix.size(); ++j) {
            genHessenberg(matrix, i, j);
        }
    }
    //outputMatrix(matrix);
}

enum qrAlgs { ORDINARY, THETHA, HESS, HESS_THETHA };


void universeQrEigenvalues(char option, vector<vector<TT>>& matrix) {
   switch (option) {
       case ORDINARY:
           for (int i = 0; i < matrix.size() - 1; ++i) {
               getNextMatrix(matrix);
           }
           getNextMatrix(matrix, true);
           break;
       case THETHA:
           for (int i = 0; i < matrix.size() - 1; ++i) {
               TT lastRow = 5.0;
               while (std::abs(lastRow) > DIVISTION_ERROR) {
                   lastRow = 0.0;
                   for(int j = matrix.size() - 1 - i; j >= 0; --j) {
                       if (matrix.size() - 1 - i != j) {
                           lastRow += matrix[matrix.size() - 1 - i][j];
                       }
                   }
                   getNextMatrix(matrix);
               }
           }
           break;
       case HESS:
           Hessenberg(matrix);
           universeQrEigenvalues(ORDINARY, matrix);
           break;
       case HESS_THETHA:
           Hessenberg(matrix);
           universeQrEigenvalues(THETHA, matrix);
           break;
       default:
           std::cerr << "wrong option // universeQrEigenvalues";
   }
}

vector<TT> qrEigenvalues(vector<vector<TT>>& matrix, vector<TT> & vec) {
    // ORDINARY, THETHA, HESS, HESS_THETHA
    universeQrEigenvalues(HESS_THETHA, matrix);

    outputOnTheScreenMatrix(matrix);

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

void getEigenvalues() {
    getMethod(&qrEigenvalues, "qrEigenvalues");
}