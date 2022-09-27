#ifndef INC_LAB1_SHARED_H
#define INC_LAB1_SHARED_H

// input file for matrix
#define IN_FILE_MATRIX "../../../test2m.txt"
// input file for vector
#define IN_FILE_VECTOR "../../../test2v.txt"
// output file for matrix
#define OUT_FILE_MATRIX "../../../outputMatrix.txt"
// output file for vector
#define OUT_FILE_VECTOR "../../../outputVector.txt"
// compare for double
#define COMPARE_RATE 10e-8
// zero division error
#define DIVISTION_ERROR 5

#endif //INC_1_SHARED_H

#include <vector>
#include <algorithm>

//  This function generates a random double in [i, j]
double GetRandomDouble(double i, double j);

void inputMatrix(std::vector<std::vector<double>>& matrix);

void outputMatrix(const std::vector<std::vector<double>>& matrix);

void outputMatrix(int amtOfVertices);

void outputVector(int amtOfElements);

void inputVector(std::vector<double>& vect);

void outputVector(const std::vector<double>& vect);

double normInfVector(const std::vector<double>& vect);

double norm1Vector(const std::vector<double>& vect);

double normInfMatrix(const std::vector<std::vector<double>>& matrix);

double norm1Matrix(const std::vector<std::vector<double>>& matrix);

std::vector<double> MultiplicationMatrixvsVector(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vect);

double normDiffer(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const std::vector<double>& x,
	double(*normVector)(const std::vector<double>&));

std::vector<std::vector<double>> transpoceMatrix(const std::vector<std::vector<double>>& matrix);

std::vector<std::vector<double>> identityMatrix(std::vector<std::vector<double>>& matrix, int size);

std::vector<std::vector<double>> inverseMatrix(std::vector<std::vector<double>>& matrix);

double condMatrix(std::vector<std::vector<double>>& A, double(*normMatrix)(const std::vector<std::vector<double>>&));

void outputOnTheScreenMatrix(const std::vector<std::vector<double>>& matrix);

void outputOnTheScreenVector(const std::vector<double>& vector);

std::vector<std::vector<double>> matrixMultiplication(const std::vector<std::vector<double>>& firstM, const std::vector<std::vector<double>>& secondM);

void disturbAndCond(std::vector<std::vector<double>>& A, std::vector<double> b, const std::vector<double>& x, double(*normVector)(const std::vector<double>&));