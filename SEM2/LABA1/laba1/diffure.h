
#pragma once

#include "main.h";
#include "matrix_operations.h"

vector<double> f(double t, vector<double> x, int dim);

vector<double> g(int dim, std::string flag, vector<double> s0, vector<double> y0, double h, const double tn);

matrix<double> Euler_explicit(vector<double> cond, int n, double h, int dim);

matrix<double> Runge_Kutta2(vector<double> cond, double h, int n, int dim);

matrix<double> Runge_Kutta4(vector<double> cond, double h, int n, int dim);

matrix<double> Adams(vector<double> cond, double h, int n, int dim);

matrix<double> prediction(vector<double> cond, double h, int n, int dim);

vector<double> Newton(std::string flag, int dim, vector<double> x0, vector<double> v, vector<double> result, double h, double tn);

matrix<double> Jacobi_matr(int dim, std::string flag, vector<double> x, vector<double> v, double h, double tn);

// Функция правой части
vector<double> f(double t, vector<double> x, int dim)
{
    vector<double> u(dim);

    // Вариант 11
  /* double r1 = 0.4; double r2 = 0.2; double b11 = 0.005;
   double b12 = 0.08; double b21 = 0.04; double b22 = 0.003;
    u[0] = r1 * x[0] - b11 * x[0] * x[0] - b12 * x[0] * x[1];
    u[1] = -r2 * x[1] - b22 * x[1] * x[1] + b21 * x[0] * x[1];*/

    // Физический маятник
    double k = 1; double m = 1;
    u[0] = x[1];
    u[1] = -k / m * x[0];
    return u;
}
//matrix<double> analit(n, vector<double>(dim))
//{
//    vector<double> u(dim);
//    u[0] = x[1];
//    u[1] = -k / m * x[0];
//    return tn * tn - tn + cos(w * tn) + sin(w * tn);// cos(w*tn);// пружина
//}

// Явный метод Эйлера 
matrix<double> Euler_explicit(vector<double> cond, int n, double h, int dim)
{
    matrix<double> y(n, vector<double>(dim));
    y[0] = cond;
    //int t = 1200;

    for (int i = 0; i < n - 1; i++)
    {
        y[i + 1] = y[i] + h * f(i * h, y[i], dim);
    }
    
    //cout << "Решение в k-м узле: " << y[t][0] << ", " << y[t][1] << endl;
    return y;
}

// Неявный метод Эйлера 
matrix<double> Euler_implicit(vector<double> cond, int n, double h, int dim)
{
    std::string flag = "Euler";
    matrix<double> y(n, vector<double>(dim));
    vector<double> temp(dim);
    int k = 0;
    y[0] = cond;
    //int t = 1200;

    for (int i = 0; i < n - 1; i++)
    {
        temp = y[i] + h * f(i * h, y[i], dim);
        y[i + 1] = Newton(flag, dim, temp, y[i], y[i + 1], h, (i + 1) * h);
    }
    //cout << "Решение в k-м узле: " << y[t][0] << ", " << y[t][1] << endl;
    return y;
}

// Симметричная схема
matrix<double> symmetric(vector<double> cond, double h, int n, int dim)
{
    std::string flag = "Symmetric";
    matrix<double> y(n, vector<double>(dim));
    vector<double> temp(dim);
    y[0] = cond;
    //int t = 1200;

    for (int i = 0; i < n - 1; i++)
    {
        temp = y[i] + h * f(i * h, y[i], dim);
        y[i + 1] = Newton(flag, dim, temp, y[i], y[i + 1], h, (i + 1) * h);
    }
    //cout << "Решение в k-м узле: " << y[t][0] << ", " << y[t][1] << endl;
    return y;
}

// Метод Рунге-Кутты 2-го порядка
matrix<double> Runge_Kutta2(vector<double> cond, double h, int n, int dim)
{
    matrix<double> y(n, vector<double>(dim));
    vector<double> k1(dim);
    vector<double> k2(dim);
    vector<double> temp(dim);
    //int t = 1200;

    y[0] = cond;
    for (int i = 0; i < n - 1; i++)
    {
        k1 = f(i * h, y[i], dim);
        temp = y[i] + h * k1;
        k2 = f((i + 1.0) * h, temp, dim);

        temp = k1 + k2;
        y[i + 1] = y[i] + h / 2.0 * temp;
    }
    //cout << "Решение в k-м узле: " << y[t][0] << ", " << y[t][1] << endl;
    return y;
}

// Метод Рунге-Кутты 4-го порядка
matrix<double> Runge_Kutta4(vector<double> cond, double h, int n, int dim)
{
    matrix<double> y(n, vector<double>(dim));
    vector<double> k1(dim);
    vector<double> k2(dim);
    vector<double> k3(dim);
    vector<double> k4(dim);
    vector<double> temp1(dim);
    vector<double> temp2(dim);
    double tau = 1.0;
    //int t = 1200;

    y[0] = cond;
    for (int i = 0; i < n - 1; i++)
    {
        // шаг = tau_{n+1}
        k1 = f(i * h, y[i], dim);
        temp1 = y[i] + 0.5 * h * k1;
        k2 = f((i + 0.5) * h, temp1, dim);
        temp1 = y[i] + 0.5 * h * k2;
        k3 = f((i + 0.5) * h, temp1, dim);
        temp1 = y[i] + 1.0 * h * k3;
        k4 = f((i + 1.0) * h, temp1, dim);

        temp1 = k1 + 2.0 * k2 + 2.0 * k3 + k4;
        temp1 = y[i] + h / 6.0 * temp1;

        // шаг = tau_{n+1/2}
        temp2 = y[i];
        for (int j = 0; j < 2; j++)
        {
            k1 = f((i + j / 2) * h, y[i], dim);
            temp2 = temp2 + tau / 4 * h * k1;
            k2 = f(((i + j / 2) + tau / 4) * h, temp2, dim);
            temp2 = temp2 + tau / 4 * h * k2;
            k3 = f(((i + j / 2) + tau / 4) * h, temp2, dim);
            temp2 = temp2 + tau / 2 * h * k3;
            k4 = f(((i + j / 2) + tau / 2) * h, temp2, dim);

            temp2 = k1 + 2.0 * k2 + 2.0 * k3 + k4;
            temp2 = y[i] + h / 6.0 * temp2;
        }
        if (norm(1 / (pow(2, 4) - 1) * (temp2 - temp1), "1") <= eps)
        {
            y[i + 1] = temp2;
        }
        else y[i + 1] = temp1;
    }
    //cout << "Решение в k-м узле: " << y[t][0] << ", " << y[t][1] << endl;
    return y;
}

// Метод Адамса 
matrix<double> Adams(vector<double> cond, double h, int n, int dim)
{
    matrix<double> y(n, vector<double>(dim));
    matrix<double> y4(n, vector<double>(dim));
    vector<double> temp(dim);
    //int t = 1200;

    y4 = Runge_Kutta4(cond, h, 4, dim);
    for (int i = 0; i < 4; i++)
    {
        y[i] = y4[i];
    }

    for (int i = 3; i < n - 1; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            temp = 55.0 * f(i * h, y[i], dim) - 59.0 * f((i - 1) * h, y[i - 1], dim);
        }
        temp = temp + 37.0 * f((i - 2) * h, y[i - 2], dim);
        temp = temp - 9.0 * f((i - 3) * h, y[i - 3], dim);
        y[i + 1] = y[i] + h / 24.0 * temp;
    }
    //cout << "Решение в k-м узле: " << y[t][0] << ", " << y[t][1] << endl;
    return y;
}

// Метод "прогноза и коррекции"
matrix<double> prediction(vector<double> cond, double h, int n, int dim)
{
    matrix<double> y(n, vector<double>(dim));
    matrix<double> y4(n, vector<double>(dim));
    vector<double> temp(dim);
    //int t = 1200;

    y4 = Runge_Kutta4(cond, h, 4, dim);
    for (int i = 0; i < 4; i++)
    {
        y[i] = y4[i];
    }

    for (int i = 3; i < n - 1; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            temp = 55.0 * f(i * h, y[i], dim) - 59.0 * f((i - 1) * h, y[i - 1], dim);
        }
        temp = temp + 37.0 * f((i - 2) * h, y[i - 2], dim);
        temp = temp - 9.0 * f((i - 3) * h, y[i - 3], dim);
        y[i + 1] = y[i] + h / 24.0 * temp;

        temp = 9.0 * f((i + 1) * h, y[i + 1], dim) + 19.0 * f(i * h, y[i], dim);
        temp = temp - 5.0 * f((i - 1) * h, y[i - 1], dim);
        temp = temp + f((i - 2) * h, y[i - 2], dim);
        y[i + 1] = y[i] + h / 24.0 * temp;
    }
    //cout << "Решение в k-м узле: " << y[t][0] << ", " << y[t][1] << endl;
    return y;
}

// Функция ..
vector<double> g(int dim, std::string flag, vector<double> s0, vector<double> y0, double h, const double tn)
{
    vector<double> ans(dim);

    if (flag == "Euler")
        for (int i = 0; i < dim; i++)
        {
            ans[i] = s0[i] - y0[i] - h * f(tn, s0, dim)[i];
        }
    if (flag == "Symmetric")
        for (int i = 0; i < dim; i++)
        {
            ans[i] = s0[i] - y0[i] - h / 2 * (f(tn, s0, dim)[i] + f(tn, y0, dim)[i]);
        }
    return ans;
}

// Метод Ньютона 
vector<double> Newton(std::string flag, int dim, vector<double> x0, vector<double> v, vector<double> result, double h, double tn)
{
    matrix<double> Jacobi(dim, vector<double>(dim));
    vector<double> ans(dim);
    vector<double> x(dim);
    vector<double> xk(dim);
    double temp; double jacobian;
    int k = 0;

    do
    {
        k++;
        x = xk;
        Jacobi = Jacobi_matr(dim, flag, x, v, h, tn);

        if (dim == 2)
        {
            jacobian = Jacobi[0][0] * Jacobi[1][1] - Jacobi[1][0] * Jacobi[0][1];
            temp = Jacobi[0][0];
            Jacobi[0][0] = Jacobi[1][1] / jacobian;
            Jacobi[1][1] = temp / jacobian;
            Jacobi[0][1] /= -jacobian;
            Jacobi[1][0] /= -jacobian;
        }
        else Jacobi = inverse(Jacobi);

        xk = Jacobi * g(dim, flag, xk, v, h, tn);
        xk = x - xk;
        x = x - xk;

        ans = xk;

    } while (norm(x, "1") > eps);
    return ans;
}

// Матрица Якоби
matrix<double> Jacobi_matr(int dim, std::string flag, vector<double> x, vector<double> v, double h, double tn)
{
    double eps1 = 1e-05;
    matrix<double> Jacobi(dim, vector<double>(dim));
    vector<double> fx = g(dim, flag, x, v, h, tn);
    vector<double> temp1(dim);
    vector<double> temp2(dim);

    for (int i = 0; i < dim; i++)
    {
        temp1 = x;
        temp1[i] += eps1;
        temp2 = g(dim, flag, temp1, v, h, tn);
        for (int j = 0; j < dim; j++)
        {
            Jacobi[j][i] = (temp2[j] - fx[j]) / eps1;
        }
    }
    return Jacobi;
}