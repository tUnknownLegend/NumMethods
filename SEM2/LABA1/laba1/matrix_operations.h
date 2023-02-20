#pragma once

#include "main.h";

template<typename T>
void read(std::string path, matrix<T>& A, int& n);

template<typename T>
void read(std::string path, vector<T>& b, int& n);

template<typename T>
void print(const matrix<T>& matr);

template<typename T>
void print(const vector<T>& vec);

template<typename T>
void identity(matrix<T>& A);

template<typename T>
matrix<T> transpose(const matrix<T>& A);

template<typename T>
matrix<T> operator+(const matrix<T>& A, const matrix<T>& B);

template<typename T>
vector<T> operator+(const vector<T>& a, const vector<T>& b);

template<typename T>
matrix<T> operator-(const matrix<T>& A, const matrix<T>& B);

template<typename T>
vector<T> operator-(const vector<T>& a, const vector<T>& b);

template<typename T>
matrix<T> operator*(const matrix<T>& A, const matrix<T>& B);

template<typename T>
vector<T> operator*(const matrix<T>& A, const vector<T>& b);

template<typename T>
matrix<T> operator*(T scalar, const matrix<T>& A);

template<typename T>
vector<T> operator*(T scalar, const vector<T>& a);

template<typename T>
matrix<T> operator-(const matrix<T>& A);

template<typename T>
T norm(const vector<T>& vec, std::string type);

template<typename T>
T norm(const matrix<T>& matr, std::string type);

template<typename T>
void diag_positive(matrix<T>& A, vector<T>& b);

template<typename T>
T norm_of_unbound(const matrix<T>& A, const vector<T>& b, const vector<T>& x, T(*norm)(const vector<T>&, std::string type));

template<typename T>
void LDU(const matrix<T>& A, matrix<T>& L, matrix<T>& D, matrix<T>& U);

void three_diag(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, double one, double two, double three);

template<typename T>
matrix<T> inverse(const matrix<T>& A);

template<typename T>
vector<T> reverse_Gauss(const matrix<T>& matr, const vector<T>& vec);

template<typename T>
vector<T> Gauss(matrix<T> matr, vector<T> vec);

template<typename T>
void read(std::string path, matrix<T>& A, int& n)
{
    std::ifstream fin;
    fin.open(path);
    double num;
    fin >> n;

    for (int i = 0; i < n; ++i)
    {
        vector<T>	vec;
        for (int j = 0; j < n; ++j)
        {
            fin >> num;
            vec.push_back(num);
        }
        A.push_back(vec);
    }
    fin.close();
}

template<typename T>
void read(std::string path, vector<T>& b, int& n)
{
    std::ifstream fin;
    fin.open(path);
    double num;
    fin >> n;

    for (int i = 0; i < n; ++i)
    {
        fin >> num;
        b.push_back(num);
    }
    fin.close();
}


template<typename T>
void print(const matrix<T>& matr)
{
    for (int i = 0; i < matr.size(); ++i)
    {
        cout << '\t';
        for (int j = 0; j < matr[0].size(); ++j)
        {
            cout << matr[i][j] << '\t';
        }
        cout << endl;
    }
}


template<typename T>
void print(const vector<T>& vec)
{
    cout << '\t';
    for (int i = 0; i < vec.size(); ++i)
    {
        cout << vec[i] << '\t';
    }
    cout << endl;
}


template<typename T>
void identity(matrix<T>& A)
{
    for (int i = 0; i < A.size(); ++i)
    {
        for (int j = 0; j < A.size(); ++j)
        {
            if (i == j) A[i][j] = 1.0;
        }
    }
}


template<typename T>
matrix<T> transpose(const matrix<T>& A)
{
    matrix<T> ans;
    for (int j = 0; j < A.size(); ++j)
    {
        vector<T> vec;
        for (int i = 0; i < A.size(); ++i)
        {
            vec.push_back(A[i][j]);
        }
        ans.push_back(vec);
    }
    return ans;
}


template<typename T>
T norm(const vector<T>& vec, std::string type)
{
    if (type == "1")
    {
        T sum = 0;
        for (int i = 0; i < vec.size(); ++i)
        {
            sum += abs(vec[i]);
        }
        return sum;
    }
    else if (type == "inf")
    {
        T max = abs(vec[0]);
        for (int i = 1; i < vec.size(); ++i)
        {
            if (max < abs(vec[i]))
            {
                max = abs(vec[i]);
            }
        }
        return max;
    }
    else
    {
        cout << "norm does noe exist";
    }
}

// Норма матрицы
template<typename T>
T norm(const matrix<T>& matr, std::string type)
{
    if (type == "1")
    {
        vector<T> vec;
        T max = 0;
        for (int j = 0; j < matr.size(); ++j)
        {
            T sum = 0;
            for (int i = 0; i < matr.size(); ++i)
            {
                sum += abs(matr[i][j]);
            }
            if (max < sum)
            {
                max = sum;
            }
        }
        return max;
    }
    else if (type == "inf")
    {
        vector<T> vec;
        T max = 0;
        for (int i = 0; i < matr.size(); ++i)
        {
            T sum = 0;
            for (int j = 0; j < matr.size(); ++j)
            {
                sum += abs(matr[i][j]);
            }
            if (max < sum)
            {
                max = sum;
            }
        }
        return max;
    }
}

// ПЕРЕГРУЗКА ОПЕРАТОРОВ

// Оператор+ для матриц
template<typename T>
matrix<T> operator+(const matrix<T>& A, const matrix<T>& B)
{
    matrix<T> ans(A.size(), vector<T>(A.size()));
    if (A.size() == B.size())
    {
        for (int i = 0; i < A.size(); ++i)
        {
            for (int j = 0; j < A.size(); ++j)
            {
                ans[i][j] = A[i][j] + B[i][j];
            }
        }
        return ans;
    }
    else
    {
        cout << "sum matrix no!!";
    }
}

// Оператор+ для векторов
template<typename T>
vector<T> operator+(const vector<T>& a, const vector<T>& b)
{
    vector<T> ans(a.size());
    if (a.size() == b.size())
    {
        for (int i = 0; i < a.size(); ++i)
        {
            ans[i] = a[i] + b[i];
        }
        return ans;
    }
    else
    {
        cout << "sum vector no";
    }
}

// Оператор- для матриц
template<typename T>
matrix<T> operator-(const matrix<T>& A, const matrix<T>& B)
{
    matrix<T> ans(A.size(), vector<T>(A.size()));
    if (A.size() == B.size())
    {
        for (int i = 0; i < A.size(); ++i)
        {
            for (int j = 0; j < A.size(); ++j)
            {
                ans[i][j] = A[i][j] - B[i][j];
            }
        }
        return ans;
    }
    else
    {
        cout << "Вычитание матриц невозможно!!!";
    }
}

// Оператор- для векторов
template<typename T>
vector<T> operator-(const vector<T>& a, const vector<T>& b)
{
    vector<T> ans(a.size());
    if (a.size() == b.size())
    {
        for (int i = 0; i < a.size(); ++i)
        {
            ans[i] = a[i] - b[i];
        }
        return ans;
    }
    else
    {
        cout << "Вычитание векторов невозможно!!!";
    }
}

// Оператор* для матриц
template<typename T>
matrix<T> operator*(const matrix<T>& A, const matrix<T>& B)
{
    matrix<T> ans(A.size(), vector<T>(A.size()));
    if (A.size() == B.size())
    {
        for (int i = 0; i < A.size(); ++i)
        {
            for (int j = 0; j < A.size(); ++j)
            {
                T s = 0;
                for (int k = 0; k < A.size(); ++k)
                {
                    s += A[i][k] * B[k][j];
                }
                ans[i][j] = s;
            }
        }
        return ans;
    }
    else
    {
        cout << "Умножение матриц невозможно!!!";
    }
}

// Оператор* для матрицы и вектора
template<typename T>
vector<T> operator*(const matrix<T>& A, const vector<T>& b)
{
    vector<T> ans(A.size());
    T s;
    if (A.size() == b.size())
    {
        for (int i = 0; i < A.size(); ++i)
        {
            s = 0;
            for (int j = 0; j < A.size(); ++j)
            {
                s += A[i][j] * b[j];
            }
            ans[i] = s;
        }
        return ans;
    }
    else
    {
        cout << "Умножение матрицы на вектор невозможно!!!";
    }
}

// Оператор* для числа и матрицы
template<typename T>
matrix<T> operator*(T scalar, const matrix<T>& A)
{
    matrix<T> ans(A.size(), vector<T>(A.size()));
    for (int i = 0; i < A.size(); ++i)
    {
        for (int j = 0; j < A.size(); ++j)
        {
            ans[i][j] = scalar * A[i][j];
        }
    }
    return ans;
}

// Оператор* для числа и вектора
template<typename T>
vector<T> operator*(T scalar, const vector<T>& a)
{
    vector<T> ans(a.size());
    for (int i = 0; i < a.size(); ++i)
    {
        ans[i] = scalar * a[i];
    }
    return ans;
}

// Оператор- (унарный) для числа и вектора
template<typename T>
matrix<T> operator-(const matrix<T>& A)
{
    matrix<T> ans(A.size(), vector<T>(A.size()));
    for (int i = 0; i < A.size(); ++i)
    {
        for (int j = 0; j < A.size(); ++j)
        {
            ans[i][j] = -A[i][j];
        }
    }
    return ans;
}

// Приведение матрицы к положительной диагонали
template<typename T>
void diag_positive(matrix<T>& A, vector<T>& b)
{
    for (int i = 0; i < A.size(); ++i)
    {
        if (A[i][i] < 0)
        {
            for (int j = 0; j < A.size(); ++j)
                A[i][j] = -A[i][j];
            b[i] = -b[i];
        }
    }
}

// Вычисление невязки 
template<typename T>
T norm_of_unbound(const matrix<T>& A, const vector<T>& b, const vector<T>& x, T(*norm)(const vector<T>&, std::string type))
{
    vector<T> razn;
    vector<T> b1 = A * x;
    for (int i = 0; i < b.size(); ++i)
    {
        razn.push_back(b[i] - b1[i]);
    }
    return norm(razn, "1");
}

// Разложение матрицы на LDU
template<typename T>
void LDU(const matrix<T>& A, matrix<T>& L, matrix<T>& D, matrix<T>& U)
{
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A.size(); j++)
        {
            if (i == j) D[i][j] = A[i][j];
            if (i > j) L[i][j] = A[i][j];
            if (i < j) U[i][j] = A[i][j];
        }
    }
}

// Инициализация трёхдиагональной матрицы
void three_diag(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, double one, double two, double three)
{
    for (int i = 0; i < a.size(); ++i)
    {
        a[i] = one;
        b[i] = two;
        c[i] = three;
        d[i] = i;
    }
    a[0] = 0; c[a.size() - 1] = 0;
}

// Вычисление обратной матрицы 
template<typename T>
matrix<T> inverse(const matrix<T>& A)
{
    matrix<T> inv;
    matrix<T> E = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1 } };
    for (int i = 0; i < A.size(); ++i)
    {
        inv.push_back(Gauss(A, E[i]));
    }
    return transpose(inv);
}

// Обратный ход метода Гаусса
template<typename T>
vector<T> reverse_Gauss(const matrix<T>& matr, const vector<T>& vec)
{
    T t;
    vector<T> x;
    for (int i = 0; i < vec.size(); ++i)
    {
        x.push_back(1);
    }

    for (int i = matr.size() - 1; i >= 0; --i)
    {
        t = vec[i];
        for (int j = i + 1; j < matr.size(); ++j)
        {
            t = t - matr[i][j] * x[j];
        }
        x[i] = t / matr[i][i];
        if (abs(x[i]) < 0.1e-10)
        {
            x[i] = 0.0;
        }
    }
    return x;
}

// Метод Гаусса
template<typename T>
vector<T> Gauss(matrix<T> matr, vector<T> vec)
{
    vector<T> ans(matr.size());
    int imax;
    T temp;
    bool flag;
    for (int k = 0; k < matr.size(); ++k)
    {
        // частичный выбор главного элемента 
        imax = k;
        for (int i = k; i < matr.size(); ++i)
        {
            if (abs(matr[i][k]) > abs(matr[imax][k]))
            {
                imax = i;
            }
        }
        flag = true;
        if (imax != k)
        {
            std::swap(matr[imax], matr[k]);
            std::swap(vec[imax], vec[k]);
        }

        // прямой ход метода Гаусса
        if (flag)
        {
            for (int i = k + 1; i < matr.size(); ++i)
            {
                temp = matr[i][k] / matr[k][k];
                vec[i] -= vec[k] * temp;
                for (int j = k; j < matr.size(); ++j)
                {
                    matr[i][j] -= matr[k][j] * temp;
                }
            }
        }
    }
    if (flag)
    {
        ans = reverse_Gauss(matr, vec);
        return ans;
    }
}
