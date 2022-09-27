#include <iostream>
#include "shared.h"

using std::vector;

vector<double> CalcGaussMethod(vector<vector<double>> matr, vector<double> vect) {

	vector<double> resultVect(vect.size(), 0.0); // инициализация вектора 

	for (int k = 0; k < vect.size(); ++k) {
		//  partial selection
		double maxValInd = k;
		for (int i = k + 1; i <vect.size(); ++i) {
			maxValInd = matr[i][k] > matr[maxValInd][k] ? i : maxValInd;
		}

		
		if (abs(matr[maxValInd][k]) > COMPARE_RATE) {
			std::swap(matr[maxValInd], matr[k]);
			std::swap(vect[maxValInd], vect[k]);

			//  straight
			for (int i = k + 1; i < vect.size(); ++i) {
				double c = (matr[i][k] / matr[k][k]);
				for (int j = k; j < vect.size(); ++j) {
					matr[i][j] -= c * matr[k][j];
				}

				vect[i] -= c * vect[k];
			}
		}
		else {
			std::cerr << "Matrix is singular";
			return {};
		}
	}

	//  reverse
	for (int i = vect.size() - 1; i >= 0; --i) {
		double sumJ = 0.0;
		for (int j = i + 1; j < vect.size(); ++j) {
			sumJ += matr[i][j] * resultVect[j];
		}

		resultVect[i] = (vect[i] - sumJ) / matr[i][i];
	}
	outputMatrix(matr);
	return resultVect;
}

vector<double> getGauss() {
	vector<vector<double>> matrix;
	vector<double> rightVect;
	inputMatrix(matrix);
	inputVector(rightVect);
	vector<double> res = CalcGaussMethod(matrix,rightVect);
	outputVector(res);
	std::cout << "Inverse matrix:" << std::endl;
	outputOnTheScreenMatrix(inverseMatrix(matrix));
	std::cout << std::endl;
	std::cout << "Matrix multiplication:" << std::endl;
	outputOnTheScreenMatrix(matrixMultiplication(matrix, inverseMatrix(matrix)));
	disturbAndCond(matrix, rightVect, res, norm1Vector);

	return res;
}