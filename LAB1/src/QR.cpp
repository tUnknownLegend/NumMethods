#include <iostream>
#include "shared.h"

using std::vector;

vector<vector<double>> matrixMultiplication(const vector<vector<double>>& firstM, const vector<vector<double>>& secondM) {
	vector<vector<double>> resMatrix(firstM.size(), vector<double>(firstM.size(), 0));

	for (int i = 0; i < firstM.size(); ++i) {
		for (int j = 0; j < secondM.size(); ++j) {
			for (int k = 0; k < firstM.size(); ++k) {
				resMatrix[i][j] += firstM[i][k] * secondM[k][j];
			}
		}
	}

	return resMatrix;
}

vector<double> MultiplicationMatrixvsVector(const vector<vector<double>>& matrix, const vector<double>& vect) {
	vector<double> resVector;
	double s;
	for (int i = 0; i < matrix.size(); ++i) {
		s = 0;
		for (int j = 0; j < matrix.size(); ++j) {
		   s += matrix[i][j] * vect[j];
		}
		resVector.push_back(s);
	}

	return resVector;
}
vector<vector<double>> transpoceMatrix(const vector<vector<double>>& matrix) {
	vector<vector<double>> resMatrix;
	vector<double> str;
	for (int j = 0; j < matrix.size(); ++j) {
		for (int i = 0; i < matrix.size(); ++i) {
			str.push_back(matrix[i][j]);
		}
		resMatrix.push_back(str);
		str.clear();
	}
	return resMatrix;
}


vector<double> CalcQRmethod() {
	vector<vector<double>> matrix;
	vector<double> rightVect;
	inputMatrix(matrix);
	inputVector(rightVect);

	vector<double> resultVect(rightVect.size(), 0.0);

	//outputMatrix(matrixMultiplication(matrix, matrix));
	//outputVector(MultiplicationMatrixvsVector(matrix, rightVect));
	outputMatrix(transpoceMatrix(matrix));
	return {};
}

vector<double> getQR() {
	vector<double> res = CalcQRmethod();
	//outputVector(res);

	return res;
}



/* vector<double> QR_decomp(const matrix<T>& A, const vector<T>& b)
{
	matrix<T> Q, R, T_vr;
	vector<T> x, add;
	T c, s;
	bool flag;
	equal(R, A);
	add = b;

	for (int i = 0; i < A.size() - 1; ++i)
	{
		for (int j = i + 1; j < A.size(); ++j)
		{
			if (abs(R[i][j]) > 0.1e-10)
			{
				c = R[i][i] / sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
				s = R[j][i] / sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
				for (int k = 0; k < A.size(); ++k)
				{
					T tm = R[i][k];
					R[i][k] = c * R[i][k] + s * R[j][k];
					R[j][k] = -s * tm + c * R[j][k];
					if (abs(R[j][k]) < 0.1e-10)
					{
						R[j][k] = 0.0;
					}
				}
				T tv = add[i];
				add[i] = c * add[i] + s * add[j];
				add[j] = -s * tv + c * add[j];
			}
		}
}*/