#include <iostream>
#include "shared.h"
#include "Gauss.h"

using std::vector;


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
	
	std::cout << "Inverse matrix:" << std::endl;
	outputOnTheScreenMatrix(inverseMatrix(matrix));
	std::cout  << std::endl;
	std::cout << "Matrix multiplication:" << std::endl;
	outputOnTheScreenMatrix(matrixMultiplication(matrix,inverseMatrix(matrix)));
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
