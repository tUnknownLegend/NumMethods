#include <iostream>
#include "shared.h"
#include "Gauss.h"

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

vector<vector<double>> identityMatrix(vector<vector<double>>& matrix, int size) {
	vector<vector<double>> resMatrix;
	vector<double> str;
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j <size; ++j) {
			if (i == j) {
				str.push_back(1.0);
			}
			else{
				str.push_back(0.0);
			}
		}
		resMatrix.push_back(str);
		str.clear();
	}
	return resMatrix;
}

  vector<vector<double>> inverseMatrix(vector<vector<double>>& matrix) {
	vector<double> res(matrix.size(),0.0);
	vector<double> str;
	vector<vector<double>> resMatrix;

	vector<vector<double>> E;
	vector<vector<double>> EE;
	E.reserve(matrix.size());
	EE =identityMatrix(E, matrix.size());

	for (int i = 0; i < matrix.size(); ++i) {
		for (int j = 0; j < matrix.size(); ++j) {
			str.push_back(EE[j][i]);
		}
		res = CalcGaussMethod(matrix,str);
		resMatrix.push_back(res);
		str.clear();
	}
	return transpoceMatrix(resMatrix);
	/*for (int i = 0; i < matrix.size(); ++i)
	{
		resMatrix.push_back(CalcGaussMethod(matrix, E[i]));
	}
	return transpoceMatrix(resMatrix);*/
}




vector<double> CalcQRmethod() {
	vector<vector<double>> matrix;
	//vector<vector<double>> Q;
	//vector<vector<double>> T_vr;
	vector<double> rightVect;
	inputMatrix(matrix);
	inputVector(rightVect);
	vector<vector<double>> R(matrix);
	double c;
	double s;
	bool flag;
	vector<double> resultVect(rightVect.size(), 0.0);
	vector<double> add(rightVect.size(), 0.0);

	for (int i = 0; i < rightVect.size() - 1; ++i) {
		for (int j = i + 1; j < rightVect.size(); ++j) {
			if (abs(R[i][j]) > COMPARE_RATE) {
				c = R[i][i] / sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
				s = R[j][i] / sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
				for (int k = 0; k < rightVect.size(); ++k) {
					double tm = R[i][k];
					R[i][k] = c * R[i][k] + s * R[j][k];
					R[j][k] = -s * tm + c * R[j][k];
				}
				double tv = add[i];
				add[i] = c * add[i] + s * add[j];
				add[j] = -s * tv + c * add[j];
			}
		}
		
	}

	vector<vector<double>> ZZ = inverseMatrix(matrix);
	outputMatrix(ZZ);
	//outputMatrix(identityMatrix(matrix));
	//outputMatrix(R);
	//outputMatrix(matrixMultiplication(matrix, matrix));
	//outputVector(MultiplicationMatrixvsVector(matrix, rightVect));
	//outputMatrix(transpoceMatrix(matrix));
	return {};
}

vector<double> getQR() {
	vector<double> res = CalcQRmethod();
	//outputVector(res);

	return res;
}
