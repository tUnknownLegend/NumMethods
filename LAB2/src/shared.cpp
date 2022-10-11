#include <fstream>
#include <iostream>
#include <random>
#include <iomanip>
#include "shared.h"

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

void inputMatrix(vector<vector<TT>> &matrix) {
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

void inputVector(vector<TT> &vect) {
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

void outputVector(const vector<TT> &vect) {
    ofstream outFile(OUT_FILE_VECTOR);
    if (!outFile.is_open()) {
        cerr << "error // output.txt open\n";
        return;
    }

    outFile << vect.size() << std::endl;

    {
        TT node = 0.0;
        for (auto &el: vect) {
            outFile << std::setprecision(8) << el << " ";
        }
        outFile << std::endl;
    }
    outFile.close();
}

void outputMatrix(const vector<vector<TT>> &matrix) {
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
        for (auto &raw: matrix) {
            for (auto &el: raw) {
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
void outputOnTheScreenMatrix(const vector<vector<TT>> &matrix) {
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix.size(); ++j) {
            cout << std::setprecision(8) << matrix[i][j] << ' ';
        }
        cout << std::endl;
    }
}

// вывод вектора на экран
void outputOnTheScreenVector(const std::vector<TT> &vector) {
    for (double i: vector) {
        cout << std::setprecision(8) << i << ' ';
    }
    cout << std::endl;
}
// Кубическая норма вектора
TT normInfVector(const vector<TT> &vect) {
    TT norm = std::abs(vect[0]);
    for (int i = 1; i < vect.size(); ++i) {
        if (norm < std::abs(vect[i]))
            norm = std::abs(vect[i]);
    }
    return norm;
}

// Октэрическая норма вектора
TT norm1Vector(const vector<TT> &vect) {
    TT norm = 0;
    for (double i: vect) {
        norm += std::abs(i);
    }
    return norm;
}

// Кубическая норма матрицы
TT normInfMatrix(const vector<vector<TT>> &matrix) {
    TT norm = 0;

    for (int j = 0; j < matrix.size(); ++j) {
        TT sum = 0;
        for (int i = 0; i < matrix.size(); ++i) {
            sum += std::abs(matrix[i][j]);
        }
        if (norm < sum)
            norm = sum;
    }
    return norm;
}

// Октаэдрическая норма матрицы
TT norm1Matrix(const vector<vector<TT>> &matrix) {
    TT norm = 0;

    for (const auto &i: matrix) {
        TT sum = 0;
        for (const auto &j: i) {
            sum += std::abs(j);
        }
        if (norm < sum)
            norm = sum;
    }
    return norm;
}

vector<TT> MultiplicationMatrixvsVector(const vector<vector<TT>> &matrix, const vector<TT> &vect) {
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
TT normDiffer(const vector<vector<TT>> &A, const vector<TT> &b, const vector<TT> &x,
              TT(*normVector)(const vector<TT> &)) {
    vector<TT> differ;
    vector<TT> b1 = MultiplicationMatrixvsVector(A, x);

    for (int i = 0; i < b.size(); ++i) {
        differ.push_back(b[i] - b1[i]);
    }
    return normVector(differ);
}

vector<vector<TT>> transpoceMatrix(const vector<vector<TT>> &matrix) {
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

vector<vector<TT>>
matrixOperations(const vector<vector<TT>> &firstM, const vector<vector<TT>> &secondM, const char &operation) {
    vector<vector<TT>> resMatrix(firstM.size(), vector<TT>(firstM.size(), 0));
    switch (operation) {
        case '*':
            for (int i = 0; i < firstM.size(); ++i) {
                for (int j = 0; j < secondM.size(); ++j) {
                    for (int k = 0; k < firstM.size(); ++k) {
                        resMatrix[i][j] += firstM[i][k] * secondM[k][j];
                    }
                }
            }
            break;
        case '+':
            for (int i = 0; i < firstM.size(); ++i) {
                for (int j = 0; j < secondM.size(); ++j) {
                    resMatrix[i][j] = firstM[i][j] + secondM[i][j];
                }
            }
            break;
        case '-':
            for (int i = 0; i < firstM.size(); ++i) {
                for (int j = 0; j < secondM.size(); ++j) {
                    resMatrix[i][j] = firstM[i][j] - secondM[i][j];
                }
            }
            break;
        default:
            std::cerr << "error";
    }
    return resMatrix;
}

void matrixDigit(const TT &digit, vector<vector<TT>> &secondM, const char &operation) {
    switch (operation) {
        case '*':
            for (auto &i: secondM) {
                for (auto &j: i) {
                    j *= digit;
                }
            }
            break;
        case '/':
            for (auto &i: secondM) {
                for (auto &j: i) {
                    j /= digit;
                }
            }
            break;
        case '+':
            for (auto &i: secondM) {
                for (auto &j: i) {
                    j += digit;
                }
            }
            break;
        case '-':
            for (auto &i: secondM) {
                for (auto &j: i) {
                    j -= digit;
                }
            }
            break;
        default:
            std::cerr << "error";
    }
}

vector<TT> matrixVectorMultiplication(const vector<vector<TT>> &firstM, const vector<TT> &secondV) {
    vector<TT> result(secondV.size(), 0);

    for (int i = 0; i < secondV.size(); ++i) {
        for (int j = 0; j < secondV.size(); ++j) {
            result[i] += firstM[i][j] * secondV[j];
        }
    }
    return result;
}

vector<TT> vectorOperation(const vector<TT> &firstV, const vector<TT> &secondV, const char &operation) {
    vector<TT> result(firstV);
    switch (operation) {
        case '*':
            for (int i = 0; i < secondV.size(); ++i) {
                result[i] *= secondV[i];
            }
            break;
        case '/':
            for (int i = 0; i < secondV.size(); ++i) {
                result[i] /= secondV[i];
            }
            break;
        case '+':
            for (int i = 0; i < secondV.size(); ++i) {
                result[i] += secondV[i];
            }
            break;
        case '-':
            for (int i = 0; i < secondV.size(); ++i) {
                result[i] -= secondV[i];
            }
            break;
        default:
            std::cerr << "error";
    }
    return result;
}

void vectorDigit(const TT &digit, vector<TT> &secondV, const char &operation) {
    switch (operation) {
        case '*':
            for (auto &i: secondV) {
                i *= digit;
            }
            break;
        case '/':
            for (auto &i: secondV) {
                i /= digit;
            }
            break;
        case '+':
            for (auto &i: secondV) {
                i += digit;
            }
            break;
        case '-':
            for (auto &i: secondV) {
                i -= digit;
            }
            break;
        default:
            std::cerr << "error";
    }
}
