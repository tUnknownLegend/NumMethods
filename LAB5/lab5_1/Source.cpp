#pragma once
#include <iostream>
#include <vector>

using namespace std;
template<typename T>
using matrix = vector<vector<T>>;

//0 <= x <= 1
template<typename T>
T test1(T x)
{
	return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}

template<typename T>
T test1Df(T x)
{
	return 0.121495 - 1.5119 * x + 5.9535 * x * x - 9.28 * x * x * x + 5. * x * x * x * x;
}


//-1 <= x <= 10
template<typename T>
T test2(T x)
{
	return sqrt(x + 1) - 1;
}

template<typename T>
T test2Df(T x)
{
	return 1. / (2. * sqrt(1. + x));
}

//0 <= x <= 1
template<typename T>
T test3(T x)
{
	return 35 * x * x * x - 67 * x * x - 3 * x + 3;
}

template<typename T>
T test3Df(T x)
{
	return -3. - 134. * x + 105. * x * x;
}






































