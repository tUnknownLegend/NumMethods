#include <fstream>
#include <iostream>
#include <random>
#include "shared.h"
#include "Gauss.h"

using std::ifstream;
using std::vector;
using std::cerr;
using std::ofstream;
using std::cout;

#define TT double

//  This function generates a random double in [i, j]
double GetRandomDouble(double i, double j) {
	std::random_device rd;  // Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(i, j);
	return dis(gen);
}

void inputMatrix(vector<vector<TT>>& matrix) {
	ifstream inFile(IN_FILE_MATRIX);
	if (!inFile.is_open()) {
		cerr << "error // input.txt open\n";
		return;
	}

	int amtOfVertices = 0;
	inFile >> amtOfVertices;
	matrix.reserve(amtOfVertices);

	{
		vector<TT> str;
		TT node = 0.0;
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

void inputVector(vector<TT>& vect) {
	ifstream inFile(IN_FILE_VECTOR);
	if (!inFile.is_open()) {
		cerr << "error // input.txt open\n";
		return;
	}

	int amtOfVertices = 0;
	inFile >> amtOfVertices;
	vect.reserve(amtOfVertices);

	{
		vector<TT> str;
		TT node = 0.0;
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
		TT node = 0.0;
		for (int j = 0; j < amtOfElements; ++j) {
			outFile << GetRandomDouble(leftBound, rightBound) << " ";
		}
		outFile << std::endl;
	}
	outFile.close();
}

void outputVector(const vector<TT>& vect) {
	ofstream outFile(OUT_FILE_VECTOR);
	if (!outFile.is_open()) {
		cerr << "error // output.txt open\n";
		return;
	}

	outFile << vect.size() << std::endl;

	{
		TT node = 0.0;
		for (auto& el : vect) {
			outFile << el << " ";
		}
		outFile << std::endl;
	}
	outFile.close();
}

void outputMatrix(const vector<vector<TT>>& matrix) {
	ofstream outFile(OUT_FILE_MATRIX);
	if (!outFile.is_open()) {
		cerr << "error // output.txt open\n";
		return;
	}

	outFile << matrix.size() << std::endl;

	{
		const int leftBound = 0;
		const int rightBound = 10;
		TT node = 0.0;
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
		TT node = 0.0;
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
void outputOnTheScreenMatrix(const vector<vector<TT>>& matrix)
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
void outputOnTheScreenVector(const std::vector<TT>& vector)
{
	for (int i = 0; i < vector.size(); ++i)
	{
		cout << vector[i] << ' ';
	}
	cout << std::endl;
}
// Кубическая норма вектора
TT normInfVector(const vector<TT>& vect)
{
	TT norm = abs(vect[0]);
	for (int i = 1; i < vect.size(); ++i)
	{
		if (norm < abs(vect[i]))
			norm = abs(vect[i]);
	}
	return norm;
}

// Октэрическая норма вектора
TT norm1Vector(const vector<TT>& vect)
{
	TT norm = 0;
	for (int i = 0; i < vect.size(); ++i)
	{
		norm += abs(vect[i]);
	}
	return norm;
}

// Кубическая норма матрицы
TT normInfMatrix(const vector<vector<TT>>& matrix)
{
	TT norm = 0;

	for (int i = 0; i < matrix.size(); ++i)
	{
		TT sum = 0;
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
TT norm1Matrix(const vector<vector<TT>>& matrix)
{
	TT norm = 0;

	for (int j = 0; j < matrix.size(); ++j)
	{
		TT sum = 0;
		for (int i = 0; i < matrix.size(); ++i)
		{
			sum += abs(matrix[i][j]);
		}
		if (norm < sum)
			norm = sum;
	}
	return norm;
}

vector<TT> MultiplicationMatrixvsVector(const vector<vector<TT>>& matrix, const vector<TT>& vect) {
	vector<TT> resVector;
	TT s;
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
TT normDiffer(const vector<vector<TT>>& A, const vector<TT>& b, const vector<TT>& x,
	 TT(*normVector)(const vector<TT>&)) {
	vector<TT> differ;
	vector<TT> b1;

	//differ.reserve(A.size());

	b1 = MultiplicationMatrixvsVector(A, x);

	for (int i = 0; i < b.size(); ++i)
	{
		differ.push_back(b[i] - b1[i]);
	}
	return normVector(differ);
}

vector<vector<TT>> transpoceMatrix(const vector<vector<TT>>& matrix) {
	vector<vector<TT>> resMatrix;
	vector<TT> str;
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
vector<vector<TT>> identityMatrix(vector<vector<TT>>& matrix, int size) {
	vector<vector<TT>> resMatrix;
	vector<TT> str;
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
vector<vector<TT>> inverseMatrix(vector<vector<TT>>& matrix) {
	vector<TT> res(matrix.size(), 0.0);
	vector<TT> str;
	vector<vector<TT>> resMatrix;

	vector<vector<TT>> E;
	vector<vector<TT>> EE;
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
vector<vector<TT>> matrixMultiplication(const vector<vector<TT>>& firstM, const vector<vector<TT>>& secondM) {
	vector<vector<TT>> resMatrix(firstM.size(), vector<TT>(firstM.size(), 0));

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
TT condMatrix(vector<vector<TT>>& A, TT(*normMatrix)(const vector<vector<TT>>&))
{
	return normMatrix(inverseMatrix(A)) * normMatrix(A);
}

// Вносим возмущение, находим решение с возмущением, сравниваем его с решением без возмущений 
void disturbAndCond(vector<vector<TT>>& A, vector<TT> b, const vector<TT>& x, TT(*normVector)(const vector<TT>&))
{
	vector<TT> db;
	for (int i = 0; i < b.size(); ++i)
	{
		db.push_back(0.01);
	}

	vector<TT> b1;
	for (int i = 0; i < b.size(); ++i)
	{
		b1.push_back(b[i] + db[i]);
	}

	vector<TT> x1(b.size(), 0.0);
	//vector<TT> x1;
	vector<TT> dx;

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

	TT deltaX = normVector(dx) / normVector(x);
	TT deltaB = normVector(db) / normVector(b);
	std::cout << "cond A >=  " << deltaX / deltaB << std::endl;
}
