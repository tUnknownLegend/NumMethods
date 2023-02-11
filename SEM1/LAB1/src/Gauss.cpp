#include <iostream>
#include "shared.h"

using std::vector;

//#define TT double


vector<TT> CalcGaussMethod(vector<vector<TT>> matr, vector<TT> vect) {

	vector<TT> resultVect(vect.size(), 0.0); // инициализация вектора 
	unsigned int mltB = 0;
	for (int k = 0; k < vect.size(); ++k) {
		//  partial selection
		TT maxValInd = k;
		for (int i = k + 1; i <vect.size(); ++i) {
			maxValInd = matr[i][k] > matr[maxValInd][k] ? i : maxValInd;
		}

		
		if (abs(matr[maxValInd][k]) > COMPARE_RATE) {
			std::swap(matr[maxValInd], matr[k]);
			std::swap(vect[maxValInd], vect[k]);

			//  straight
			for (int i = k + 1; i < vect.size(); ++i) {
				TT c = (matr[i][k] / matr[k][k]);
				for (int j = k; j < vect.size(); ++j) {
					++mltB;
					matr[i][j] -= c * matr[k][j];
				}

				++mltB;
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
		TT sumJ = 0.0;
		for (int j = i + 1; j < vect.size(); ++j) {
			++mltB;
			sumJ += matr[i][j] * resultVect[j];
		}

		resultVect[i] = (vect[i] - sumJ) / matr[i][i];
	}
	outputMatrix(matr);

	//std::cout << "\nAmount of mltB: " << mltB << "\n";
	return resultVect;
}

vector<TT> getGauss() {
	vector<vector<TT>> matrix;
	vector<TT> rightVect;
	inputMatrix(matrix);
	inputVector(rightVect);
	vector<TT> res = CalcGaussMethod(matrix,rightVect);
	std::cout << "Gauss method:\n";
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

	for (TT i = 0.001; i <= 0.1; i*=3) {
		std::cout << "i = " << i << "\n";
		disturbAndCond(matrix, rightVect, res, norm1Vector, i);
	}

	return res;
}