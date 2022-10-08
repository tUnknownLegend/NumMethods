#ifndef INC_LAB1_SHARED_H
#define INC_LAB1_SHARED_H

// input file for matrix
#define IN_FILE_MATRIX "../../../test5m.txt"
// input file for vector
#define IN_FILE_VECTOR "../../../test5v.txt"
// output file for matrix
#define OUT_FILE_MATRIX "../../../outputMatrix.txt"
// output file for vector
#define OUT_FILE_VECTOR "../../../outputVector.txt"
// compare for double
#define COMPARE_RATE 10e-8
// zero division error
#define DIVISTION_ERROR 5

#include <vector>
#include <algorithm>

#define TT double

//  This function generates a random TT in [i, j]
double GetRandomDouble(double i, double j);

void inputMatrix(std::vector<std::vector<TT>>& matrix);

void outputMatrix(const std::vector<std::vector<TT>>& matrix);

void outputMatrix(int amtOfVertices);

void outputVector(int amtOfElements);

void inputVector(std::vector<TT>& vect);

void outputVector(const std::vector<TT>& vect);

TT normInfVector(const std::vector<TT>& vect);

TT norm1Vector(const std::vector<TT>& vect);

TT normInfMatrix(const std::vector<std::vector<TT>>& matrix);

TT norm1Matrix(const std::vector<std::vector<TT>>& matrix);

std::vector<TT> MultiplicationMatrixvsVector(const std::vector<std::vector<TT>>& matrix, const std::vector<TT>& vect);

TT normDiffer(const std::vector<std::vector<TT>>& A, const std::vector<TT>& b, const std::vector<TT>& x,
	TT(*normVector)(const std::vector<TT>&));

std::vector<std::vector<TT>> transpoceMatrix(const std::vector<std::vector<TT>>& matrix);

std::vector<std::vector<TT>> identityMatrix(int size);

std::vector<std::vector<TT>> inverseMatrix(std::vector<std::vector<TT>>& matrix);

TT condMatrix(std::vector<std::vector<TT>>& A, TT(*normMatrix)(const std::vector<std::vector<TT>>&));

void outputOnTheScreenMatrix(const std::vector<std::vector<TT>>& matrix);

void outputOnTheScreenVector(const std::vector<TT>& vector);

std::vector<std::vector<TT>> matrixMultiplication(const std::vector<std::vector<TT>>& firstM, const std::vector<std::vector<TT>>& secondM);

void disturbAndCond(std::vector<std::vector<TT>>& A, std::vector<TT> b, const std::vector<TT>& x, TT(*normVector)(const std::vector<TT>&), const TT diffVal);

#endif //INC_1_SHARED_H
