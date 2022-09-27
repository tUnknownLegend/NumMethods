#include <iostream>
#include "shared.h"
#include "Gauss.h"

using std::vector;

#define TT double

vector<TT> CalcQRmethod() {
	vector<vector<TT>> matrix;
	vector<vector<TT>> Q;
	vector<TT> rightVect;
	inputMatrix(matrix);
	inputVector(rightVect);
	vector<vector<TT>> R(matrix);
	vector<vector<TT>> T_vr;
	T_vr.reserve(matrix.size());
	vector<vector<TT>> T = identityMatrix(T_vr, matrix.size());
	TT c;
	TT s;
	bool flag = true;
	vector<TT> resultVect(rightVect.size(), 0.0);
	vector<TT> add(rightVect);

	for (int i = 0; i < rightVect.size() - 1; ++i) {
		for (int j = i + 1; j < rightVect.size(); ++j) {
			if (abs(R[i][j]) > COMPARE_RATE) {
				c = R[i][i] / sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
				s = R[j][i] / sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
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
		for (int j = i + 1; j < rightVect.size(); ++j)
		{
			sumJ += R[i][j] * resultVect[j];
		}
		resultVect[i] = (add[i] - sumJ) / R[i][i];

	}

	for (int i = 0; i < matrix.size(); ++i)
	{
		if (std::abs(R[i][i]) < COMPARE_RATE)
		{
			flag = 0;
			break;
		}
		else
			flag = 1;
	}

	if (flag)
	{
		Q = transpoceMatrix(T);

		std::cout << "Q:" << std::endl;
		outputOnTheScreenMatrix(Q);
		std::cout << "R:" << std::endl;
		outputOnTheScreenMatrix(R);

	}
	else {
			std::cerr << "Matrix is singular";
			return {};
		}

	//outputMatrix(identityMatrix(matrix));
	//outputMatrix(R);
	//outputMatrix(matrixMultiplication(matrix, matrix));
	//outputVector(MultiplicationMatrixvsVector(matrix, rightVect));
	//outputMatrix(transpoceMatrix(matrix));
	return {resultVect};
}

vector<TT> getQR() {
	vector<vector<TT>> matrix;
	vector<TT> rightVect;
	inputMatrix(matrix);
	inputVector(rightVect);
	vector<TT> res = CalcQRmethod();
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
	disturbAndCond(matrix, rightVect, res, norm1Vector);
	return res;
}
