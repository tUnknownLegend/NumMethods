#include <fstream>
#include <iostream>
#include <random>
#include <iomanip>
#include "QR.h"
#include "shared.h"
#include "Gauss.h"

using std::ifstream;
using std::vector;
using std::cerr;
using std::ofstream;
using std::cout;

//#define TT double

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
			outFile << std::setprecision(8) << GetRandomDouble(leftBound, rightBound) << " ";
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
			outFile << std::setprecision(8) << el << " ";
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
				outFile << std::setprecision(8) << el << " ";
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
				outFile << std::setprecision(8) << GetRandomDouble(leftBound, rightBound) << " ";
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
			cout << std::setprecision(8) << matrix[i][j] << ' ';
		}
		cout << std::endl;
	}
}

// вывод вектора на экран
void outputOnTheScreenVector(const std::vector<TT>& vector)
{
	for (int i = 0; i < vector.size(); ++i)
	{
		cout << std::setprecision(8) << vector[i] << ' ';
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
	vector<TT> b1 = MultiplicationMatrixvsVector(A, x);

	for (int i = 0; i < b.size(); ++i) {
		differ.push_back(b[i] - b1[i]);
	}
	return normVector(differ);
}

vector<vector<TT>> transpoceMatrix(const vector<vector<TT>>& matrix) {
	vector<vector<TT>> resMatrix(matrix);
	for (int j = 0; j < matrix.size(); ++j) {
		for (int i = 0; i < matrix.size(); ++i) {
			resMatrix[j][i] = matrix[i][j];
		}
	}
	return resMatrix;
}

// Единичная матрица
vector<vector<TT>> identityMatrix(int size) {
	vector<vector<TT>> resMatrix(size, vector<TT>(size, 0.0));
	for (int i = 0; i < size; ++i) {
		resMatrix[i][i] = 1.0;
	}
	return resMatrix;
}

void CalcR(const vector<vector<TT>>& matrix, const vector<TT> rightVect, vector<vector<TT>>& R) {
	TT c;
	TT s;
	unsigned int coounter = 0;
	unsigned int mltB = 0;

	for (int i = 0; i < rightVect.size(); ++i) {  // - 1
		for (int j = i + 1; j < rightVect.size(); ++j) {
			if (abs(R[i][j]) > COMPARE_RATE) {
				coounter += 6;
				TT tempres = sqrt(pow(R[i][i], 2) + pow(R[j][i], 2));
				c = R[i][i] / tempres;
				s = R[j][i] / tempres;
				for (int k = 0; k < rightVect.size(); ++k) {
					TT tr = R[i][k];
					coounter += 6;
					R[i][k] = c * R[i][k] + s * R[j][k];
					R[j][k] = -s * tr + c * R[j][k];

				}
			}
		}
	}

	for (int i = 0; i < matrix.size(); ++i) {
		if (std::abs(R[i][i]) < COMPARE_RATE) {
			std::cerr << "Matrix is singular";
			return;
		}
	}

	return;
}

vector<TT> CalcPartQR(vector<vector<TT>>& matrix, vector<TT> rightVect, const vector<vector<TT>>& R) {
	TT c;
	TT s;
	vector<TT> resultVect(rightVect.size(), 0.0);
	vector<TT> add(rightVect);
	unsigned int coounter = 0;
	unsigned int mltB = 0;

	for (int i = 0; i < rightVect.size(); ++i) {  // - 1
		for (int j = i + 1; j < rightVect.size(); ++j) {
			if (abs(R[i][j]) > COMPARE_RATE) {
				coounter += 6;
				TT tempres = sqrt(pow(R[i][i], 2) + pow(R[j][i], 2));
				c = R[i][i] / tempres;
				s = R[j][i] / tempres;

				TT tv = add[i];
				mltB += 4;
				add[i] = c * add[i] + s * add[j];
				add[j] = -s * tv + c * add[j];
			}
		}
	}

	for (int i = rightVect.size() - 1; i >= 0; --i) {
		TT sumJ = 0.0;
		for (int j = i + 1; j < rightVect.size(); ++j) {
			++coounter;
			sumJ += R[i][j] * resultVect[j];
		}
		coounter += 2;
		resultVect[i] = (add[i] - sumJ) / R[i][i];
	}

	resultVect = CalcGaussMethod(matrix, rightVect);
	for (int i = 0; i < matrix.size(); ++i) {
		if (std::abs(R[i][i]) < COMPARE_RATE) {
			std::cerr << "Matrix is singular";
			return {};
		}
	}

	return resultVect;
}


// Обратная матрица
vector<vector<TT>> inverseMatrix(vector<vector<TT>> &matrix) {
    vector<TT> res(matrix.size(), 0.0);
    vector<TT> str;
    vector<vector<TT>> resMatrix;
    vector<vector<TT>> EE;
    EE = identityMatrix(matrix.size());

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
void disturbAndCond(vector<vector<TT>>& A, vector<TT> b, const vector<TT>& x, TT(*normVector)(const vector<TT>&), const TT diffVal)
{
	vector<TT> db (A.size(), 0.0);
	db[0] = diffVal;

	for (int i = 0; i < b.size(); ++i)
	{
		db.push_back(diffVal);
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
