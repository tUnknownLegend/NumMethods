#include <fstream>
#include <iostream>
#include <random>
#include "shared.h"

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
<<<<<<< Updated upstream
}
=======
}

// вывод матрицы на экран
void outputOnTheScreenMatrix(const vector<vector<double>>& matrix) {
	for (int i = 0; i < matrix.size(); ++i) {
		for (int j = 0; j < matrix.size(); ++j) {
			cout << matrix[i][j] << ' ';
		}
		cout << std::endl;
	}
}

// вывод вектора на экран
void outputOnTheScreenVector(const std::vector<double>& vector) {
	for (int i = 0; i < vector.size(); ++i) {
		cout << vector[i] << ' ';
	}
	cout << std::endl;
}

// Кубическая норма вектора
double normInfVector(const vector<double>& vect) {
	double norm = abs(vect[0]);
	for (int i = 1; i < vect.size(); ++i) {
		if (norm < abs(vect[i]))
			norm = abs(vect[i]);
	}
	return norm;
}

// Октэрическая норма вектора
double norm1Vector(const vector<double>& vect) {
	double norm = 0;
	for (int i = 0; i < vect.size(); ++i) {
		norm += abs(vect[i]);
	}
	return norm;
}

// Кубическая норма матрицы
double normInfMatrix(const vector<vector<double>>& matrix) {
	double norm = 0;

	for (int i = 0; i < matrix.size(); ++i) {
		double sum = 0;
		for (int j = 0; j < matrix.size(); ++j) {
			sum += abs(matrix[i][j]);
		}
		if (norm < sum)
			norm = sum;
	}
	return norm;
}

// Октаэдрическая норма матрицы
double norm1Matrix(const vector<vector<double>>& matrix) {
	double norm = 0;

	for (int j = 0; j < matrix.size(); ++j) {
		double sum = 0;
		for (int i = 0; i < matrix.size(); ++i) {
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
	vector<double> b1 = MultiplicationMatrixvsVector(A, x);;

	for (int i = 0; i < b.size(); ++i) {
		differ.push_back(b[i] - b1[i]);
	}
	return normVector(differ);
}

vector<vector<double>> transpoceMatrix(const vector<vector<double>>& matrix) {
	vector<vector<double>> resMatrix(size, std::vector<int>(size, 0.0));
	for (int j = 0; j < matrix.size(); ++j) {
		for (int i = 0; i < matrix.size(); ++i) {
			resMatrix[j][i] = matrix[i][j]
		}
	}
	return resMatrix;
}

// Единичная матрица
vector<vector<double>> identityMatrix(const vector<vector<double>>& matrix, int size) {
	vector<vector<double>> resMatrix(size, std::vector<int>(size, 0.0));
	for (int i = 0; i < size; ++i) {
		resMatrix[i][i] = 1.0;
	}
	return resMatrix;
}

// Обратная матрица
vector<vector<double>> inverseMatrix(const vector<vector<double>>& matrix) {
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
// Вычисление числа обусловленности
double condMatrix(vector<vector<double>>& A, double(*normMatrix)(const vector<vector<double>>&))
{
	return normMatrix(inverseMatrix(A)) * normMatrix(A);
}


>>>>>>> Stashed changes
