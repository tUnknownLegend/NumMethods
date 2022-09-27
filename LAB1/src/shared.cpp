﻿#include <fstream>
#include <iostream>
#include <random>
#include "shared.h"
#include "Gauss.h"

using std::ifstream;
using std::vector;
using std::cerr;
using std::ofstream;
using std::cout;

//  This function generates a random double in [i, j]
double GetRandomDouble(double i, double j) {
	std::random_device rd;  // Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(i, j);
	return dis(gen);
}

void inputMatrix(vector<vector<double>>& matrix) {
	ifstream inFile(IN_FILE_MATRIX);
	if (!inFile.is_open()) {
		cerr << "error // input.txt open\n";
		return;
	}

	int amtOfVertices = 0;
	inFile >> amtOfVertices;
	matrix.reserve(amtOfVertices);

	{
		vector<double> str;
		double node = 0.0;
		for (int i = 0; i < amtOfVertices; ++i) {

			for (int j = 0; j < amtOfVertices; ++j) {
				inFile >> node;
				str.push_back(node);
			}
			matrix.push_back(str);
			str.clear();
		}
	}
	inFile.close();
}

void inputVector(vector<double>& vect) {
	ifstream inFile(IN_FILE_VECTOR);
	if (!inFile.is_open()) {
		cerr << "error // input.txt open\n";
		return;
	}

	int amtOfVertices = 0;
	inFile >> amtOfVertices;
	vect.reserve(amtOfVertices);

	{
		vector<double> str;
		double node = 0.0;
		for (int j = 0; j < amtOfVertices; ++j) {
			inFile >> node;
			str.push_back(node);
		}
		vect = std::move(str);
	}
	inFile.close();
}

void outputVector(int amtOfElements) {
	ofstream outFile(OUT_FILE_VECTOR);
	if (!outFile.is_open()) {
		cerr << "error // output.txt open\n";
		return;
	}

	outFile << amtOfElements << std::endl;

	{
		const int leftBound = 1;
		const int rightBound = 10;
		double node = 0.0;
		for (int j = 0; j < amtOfElements; ++j) {
			outFile << GetRandomDouble(leftBound, rightBound) << " ";
		}
		outFile << std::endl;
	}
	outFile.close();
}

void outputVector(const vector<double>& vect) {
	ofstream outFile(OUT_FILE_VECTOR);
	if (!outFile.is_open()) {
		cerr << "error // output.txt open\n";
		return;
	}

	outFile << vect.size() << std::endl;

	{
		double node = 0.0;
		for (auto& el : vect) {
			outFile << el << " ";
		}
		outFile << std::endl;
	}
	outFile.close();
}

void outputMatrix(const vector<vector<double>>& matrix) {
	ofstream outFile(OUT_FILE_MATRIX);
	if (!outFile.is_open()) {
		cerr << "error // output.txt open\n";
		return;
	}

	outFile << matrix.size() << std::endl;

	{
		const int leftBound = 0;
		const int rightBound = 10;
		double node = 0.0;
		for (auto& raw : matrix) {
			for (auto& el : raw) {
				outFile << el << " ";
			}
			outFile << std::endl;
		}
	}
	outFile.close();
}

void outputMatrix(int amtOfVertices) {
	ofstream outFile(OUT_FILE_MATRIX);
	if (!outFile.is_open()) {
		cerr << "error // output.txt open\n";
		return;
	}

	outFile << amtOfVertices << std::endl;

	{
		const int leftBound = 0;
		const int rightBound = 10;
		double node = 0.0;
		for (int i = 0; i < amtOfVertices; ++i) {

			for (int j = 0; j < amtOfVertices; ++j) {
				outFile << GetRandomDouble(leftBound, rightBound) << " ";
			}
			outFile << std::endl;
		}
	}
	outFile.close();
}

// вывод матрицы на экран
void outputOnTheScreenMatrix(const vector<vector<double>>& matrix)
{
	for (int i = 0; i < matrix.size(); ++i)
	{
		for (int j = 0; j < matrix.size(); ++j)
		{
			cout << matrix[i][j] << ' ';
		}
		cout << std::endl;
	}
}

// вывод вектора на экран
void outputOnTheScreenVector(const std::vector<double>& vector)
{
	for (int i = 0; i < vector.size(); ++i)
	{
		cout << vector[i] << ' ';
	}
	cout << std::endl;
}
// Кубическая норма вектора
double normInfVector(const vector<double>& vect)
{
	double norm = abs(vect[0]);
	for (int i = 1; i < vect.size(); ++i)
	{
		if (norm < abs(vect[i]))
			norm = abs(vect[i]);
	}
	return norm;
}

// Октэрическая норма вектора
double norm1Vector(const vector<double>& vect)
{
	double norm = 0;
	for (int i = 0; i < vect.size(); ++i)
	{
		norm += abs(vect[i]);
	}
	return norm;
}

// Кубическая норма матрицы
double normInfMatrix(const vector<vector<double>>& matrix)
{
	double norm = 0;

	for (int i = 0; i < matrix.size(); ++i)
	{
		double sum = 0;
		for (int j = 0; j < matrix.size(); ++j)
		{
			sum += abs(matrix[i][j]);
		}
		if (norm < sum)
			norm = sum;
	}
	return norm;
}

// Октаэдрическая норма матрицы
double norm1Matrix(const vector<vector<double>>& matrix)
{
	double norm = 0;

	for (int j = 0; j < matrix.size(); ++j)
	{
		double sum = 0;
		for (int i = 0; i < matrix.size(); ++i)
		{
			sum += abs(matrix[i][j]);
		}
		if (norm < sum)
			norm = sum;
	}
	return norm;
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

// Норма невязки 
double normDiffer(const vector<vector<double>>& A, const vector<double>& b, const vector<double>& x,
	 double(*normVector)(const vector<double>&)) {
	vector<double> differ;
	vector<double> b1;

	//differ.reserve(A.size());

	b1 = MultiplicationMatrixvsVector(A, x);

	for (int i = 0; i < b.size(); ++i)
	{
		differ.push_back(b[i] - b1[i]);
	}
	return normVector(differ);
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

// Единичная матрица
vector<vector<double>> identityMatrix(vector<vector<double>>& matrix, int size) {
	vector<vector<double>> resMatrix;
	vector<double> str;
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			if (i == j) {
				str.push_back(1.0);
			}
			else {
				str.push_back(0.0);
			}
		}
		resMatrix.push_back(str);
		str.clear();
	}
	return resMatrix;
}

// Обратная матрица
vector<vector<double>> inverseMatrix(vector<vector<double>>& matrix) {
	vector<double> res(matrix.size(), 0.0);
	vector<double> str;
	vector<vector<double>> resMatrix;

	vector<vector<double>> E;
	vector<vector<double>> EE;
	E.reserve(matrix.size());
	EE = identityMatrix(E, matrix.size());

	for (int i = 0; i < matrix.size(); ++i) {
		for (int j = 0; j < matrix.size(); ++j) {
			str.push_back(EE[j][i]);
		}
		res = CalcGaussMethod(matrix, str);
		resMatrix.push_back(res);
		str.clear();
	}
	return transpoceMatrix(resMatrix);
}

//Умножение матриц
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

// Вычисление числа обусловленности
double condMatrix(vector<vector<double>>& A, double(*normMatrix)(const vector<vector<double>>&))
{
	return normMatrix(inverseMatrix(A)) * normMatrix(A);
}

// Вносим возмущение, находим решение с возмущением, сравниваем его с решением без возмущений 
void disturbAndCond(vector<vector<double>>& A, vector<double> b, const vector<double>& x, double(*normVector)(const vector<double>&))
{
	vector<double> db;
	for (int i = 0; i < b.size(); ++i)
	{
		db.push_back(0.01);
	}

	vector<double> b1;
	for (int i = 0; i < b.size(); ++i)
	{
		b1.push_back(b[i] + db[i]);
	}

	vector<double> x1(b.size(), 0.0);
	//vector<double> x1;
	vector<double> dx;

	std::cout << "Vector b with disturbance: ";
	outputOnTheScreenVector(b1);
	std::cout << std::endl;

	x1 = CalcGaussMethod(A, b1);
	std::cout << "Vector x with disturbance: ";
	outputOnTheScreenVector(x1);
	std::cout << std::endl;

	//cout << "\nНорма вектора невязки (с возмущением): " << normResidual(A, b1, x1, norm_1_vec);

	std::cout << "Comparison of solutions: ";
	for (int i = 0; i < x.size(); ++i)
	{
		dx.push_back(x1[i] - x[i]);
	}
	outputOnTheScreenVector(dx);
	std::cout << std::endl;

	double deltaX = normVector(dx) / normVector(x);
	double deltaB = normVector(db) / normVector(b);
	std::cout << "cond A >=  " << deltaX / deltaB << std::endl;
}
