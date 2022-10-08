#include <iostream>
#include "shared.h"
#include "Gauss.h"

using std::vector;

//#define TT double

vector<TT> CalcQRmethod(vector<vector<TT>>& matrix, vector<TT> rightVect) {
	vector<vector<TT>> R(matrix);
	vector<vector<TT>> T = identityMatrix(matrix.size());
	TT c;
	TT s;
	vector<TT> resultVect(rightVect.size(), 0.0);
	vector<TT> add(rightVect);
	unsigned int coounter = 0;
	unsigned int mltB = 0;

	for (int i = 0; i < rightVect.size(); ++i) {  // - 1
		for (int j = i + 1; j < rightVect.size(); ++j) {
			if (abs(R[i][j]) > COMPARE_RATE) {
				coounter += 6;
				TT tempres = sqrt(pow(R[i][i], 2) + pow(R[j][i], 2));
				c = R[i][i] / tempres;
				s = R[j][i] / tempres;
				for (int k = 0; k < rightVect.size(); ++k) {
					TT tr = R[i][k];
					coounter += 6;
					R[i][k] = c * R[i][k] + s * R[j][k];
					R[j][k] = -s * tr + c * R[j][k];

					coounter += 6;
					TT tq = T[i][k];
					T[i][k] = c * T[i][k] + s * T[j][k];
					T[j][k] = -s * tq + c * T[j][k];
				}
				coounter += 6;
				TT tv = add[i];
				mltB += 4;
				add[i] = c * add[i] + s * add[j];
				add[j] = -s * tv + c * add[j];
			}
		}	
	}

	for (int i = rightVect.size() - 1; i >= 0; --i) {
		TT sumJ = 0.0;
		for (int j = i + 1; j < rightVect.size(); ++j) {
			++coounter;
			sumJ += R[i][j] * resultVect[j];
		}
		coounter += 2;
		resultVect[i] = (add[i] - sumJ) / R[i][i];
	}

	for (int i = 0; i < matrix.size(); ++i) {
		if (std::abs(R[i][i]) < COMPARE_RATE) {
			std::cerr << "Matrix is singular";
			return {};
		}
	}

	vector<vector<TT>> Q = transpoceMatrix(T);
	std::cout << "Q:" << std::endl;
	outputOnTheScreenMatrix(Q);
	std::cout << "R:" << std::endl;
	outputOnTheScreenMatrix(R);

	//std::cout << "\nAmount of mult: " << coounter << "\n";
	//std::cout << "\nAmount of mltB: " << mltB << "\n";

	return resultVect;
}

vector<TT> getQR() {
	vector<vector<TT>> matrix;
	vector<TT> rightVect;
	inputMatrix(matrix);
	inputVector(rightVect);
	vector<TT> res = CalcQRmethod(matrix, rightVect);
	std::cout << "\n------------------------------------------\nQR method:\n";
	outputVector(res);
	std::cout << "Inverse matrix:" << std::endl;
	outputOnTheScreenMatrix(inverseMatrix(matrix));
	std::cout << std::endl;
	std::cout << "Matrix multiplication:" << std::endl;
	outputOnTheScreenMatrix(matrixMultiplication(matrix, inverseMatrix(matrix)));
	std::cout << std::endl;
	std::cout << "Residual norm : " << normDiffer(matrix, rightVect, res, norm1Vector) << std::endl;
	std::cout << std::endl;
	std::cout << "condA_1 = " << condMatrix(matrix, norm1Matrix) << std::endl;
	std::cout << std::endl;
	std::cout << "condA_inf = " << condMatrix(matrix, normInfMatrix) << std::endl;
	std::cout << std::endl;
	
	for (TT i = 0.001; i <= 0.1; i *= 3) {
		std::cout << "i = " << i << "\n";
		disturbAndCond(matrix, rightVect, res, norm1Vector, i);
	}

	return res;
}
