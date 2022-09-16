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
	//vector<vector<double>> R; 
	//vector<vector<double>> Q;
	//vector<vector<double>> T_vr;
	vector<double> rightVect;
	inputMatrix(matrix);
	inputVector(rightVect);
	double c;
	double s;
	bool flag;
	vector<double> resultVect(rightVect.size(), 0.0);
	vector<double> add(rightVect.size(), 0.0);

	for (int i = 0; i < rightVect.size() - 1; ++i) {
		for (int j = i + 1; j < rightVect.size(); ++j) {
			if (abs(matrix[i][j]) > COMPARE_RATE) {
				c = matrix[i][i] / sqrt(matrix[i][i] * matrix[i][i] + matrix[j][i] * matrix[j][i]);
				s = matrix[j][i] / sqrt(matrix[i][i] * matrix[i][i] + matrix[j][i] * matrix[j][i]);
				for (int k = 0; k < rightVect.size(); ++k) {
					double tm = matrix[i][k];
					matrix[i][k] = c * matrix[i][k] + s * matrix[j][k];
					matrix[j][k] = -s * tm + c * matrix[j][k];
				}
				double tv = add[i];
				add[i] = c * add[i] + s * add[j];
				add[j] = -s * tv + c * add[j];
			}
		}
		
	}
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
