#ifndef LAB5_TEST_H
#define LAB5_TEST_H

#include <array>
#include <vector>
#include <algorithm>

double phi1(double x);

double phi2(double x);

double phi3(double x);

double phi4(double x);

double phi5(double x);

double phi6(double x);

double phi7(double x);

double phi8(double x);

double phi9(double x);

double phi10(double x);

double phi11(double x);

double psi1(double s);

double psi2(double s);

double psi3(double s);

double psi4(double s);

double psi5(double s);

double psi6(double s);

double psi7(double s);

double psi8(double s);

double psi9(double s);

double psi10(double s);

double psi11(double s);

// Массив функций вырожденного ядра
std::vector<std::function<double(double)>> phi = {phi1, phi2, phi3,
                                                  phi4, phi5, phi6,
                                                  phi7, phi8, phi9,
                                                  phi10, phi11};

// Массив функций вырожденного ядра
std::vector<std::function<double(double)>> psi = {psi1, psi2, psi3,
                                                  psi4, psi5, psi6,
                                                  psi7, psi8, psi9,
                                                  psi10, psi11};

// Правая часть
double f(double x);

// Ядро оператора Фредгольма
double K(double x, double s);

// Область интегрирования
const std::array<double, 2> area = {0.0, 1.0};

// Коэффициент перед интегралом
const double lambda = 0.5;

// Правая часть для задачи сингулярного ядра
double F(double phi);

// Сингулярное ядро
std::array<double, 2> Q(std::array<double, 2> r, std::array<double, 2> p);

#endif //LAB5_TEST_H
