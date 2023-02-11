#include <algorithm>
#include "Source.cpp"
#include "Source1.cpp"
#include "shared.h"

// локализация корней
template <typename T, typename F>
vector<pair<T, T>> localization(F fun, T a, T b, size_t n)
{
	T h = (b - a) / n;
	vector<T> setka;
	for (size_t i = 0; i < n; i++) {
		setka.push_back(a + h * i);
	}
	setka.push_back(b);
	pair<T, T> local;
	vector<pair<T, T>> locals;
	vector<T> roots;
	size_t i = 0;
	while (i < n) {
		if (fun(setka[i]) * fun(setka[i + 1]) < 0.0) {
			local.first = setka[i];
			local.second = setka[i + 1];
			locals.push_back(local);
		}
		++i;
	}

	return locals;
}

// Метод бисекции
template<typename T, typename F>
T bisection(F fun, T a, T b, T eps, size_t& iterCount) {
	T ak = a;
	T bk = b;

	T fak = fun(ak);
	T fbk = fun(bk);
	T xk = (a + b) / 2;
	vector<T> X;
	X.push_back(xk);
	while (abs(bk - ak) > 2 * eps) {
		T fxk = fun(xk);
		if (fak * fxk <= 0) {
			bk = xk;
			fbk = fxk;
		}
		else {
			ak = xk;
			fak = fxk;
		}

		xk = (ak + bk) / 2;
		X.push_back(xk);
		++iterCount;
	}
	return xk;
}

template<typename T, typename F>
vector<T> bisection(F fun, vector<pair<T, T>> intervals, T eps) {
	vector<T> res;
	res.reserve(intervals.size());
	size_t iterCount = 0;
	for (pair<T, T>& interval : intervals) {
		T a = interval.first;
		T b = interval.second;
		T x = bisection(fun, a, b, eps, iterCount);
		res.push_back(x);
	}
	cout << "iterations:" << iterCount << endl;
	return res;
}

// Классический метод Ньютона 
// diff - производная функции fun (передаем в функцию либо аналитически вычисленную производную, либо численно вычисленную производную)
template<typename T, typename F, typename Diff>
T Newton(F fun, Diff diff, T a, T b, T x0, T eps, size_t& iterCount, size_t p = 1) {
	T xk = x0;
	T coef;
	vector<T> X;
	X.push_back(xk);
	do {
		T newXk = xk - p * fun(xk) / diff(xk);
		coef = abs(newXk - xk);
		xk = newXk;
		X.push_back(xk);
		++iterCount;
	} while (coef > eps);
	degreeX(X, 0.0);
	return xk;
}



// Модификация метода Ньютона с защитой от выхода за границу
template<typename T, typename F, typename Diff>
T NewtonModification(F fun, Diff diff, T a, T b, T x0, T eps, size_t& iterCount) {
	T fa = fun(a);
	T fb = fun(b);
	T xk = x0;
	T coef;
	//size_t iterCount = 0;

	do {
		T newXk;
		if ((iterCount + 1) % 7!= 0)
			newXk = xk - fun(xk) / diff(xk);
		else {
			// Защита от зацикливания (просто раз в 7 итераций применяем
			// другой метод)
			T fxk = fun(xk);


			// метод бисекции
			/*if (fa * fxk <= 0)
				newXk = (a + xk) / 2;
			else
				newXk = (xk + b) / 2;*/

			// модифицированный метод хорд
			if (fa * fxk <= 0)
				newXk = xk - fxk * (a - xk) / (fa - fxk);
			else
				newXk = xk - fxk * (b - xk) / (fb - fxk);
		}

		// Защита модифицированным методом бисекции
		if (newXk < a || newXk > b) {
			T fxk = fun(xk);

			// метод бисекции
			/*if (fa * fxk <= 0)
				newXk = (a + xk) / 2;
			else
				newXk = (xk + b) / 2;
				*/

				// модифицированный метод хорд
			if (fa * fxk <= 0)
				newXk = xk - fxk * (a - xk) / (fa - fxk);
			else
				newXk = xk - fxk * (b - xk) / (fb - fxk);

		}

		coef = abs(newXk - xk);
		xk = newXk;
		iterCount++;
	} while (coef > eps);

	cout << "iterations: " << iterCount << endl;

	return xk;
}




// Метод Ньютона
// diff - производная функции fun
template<typename T, typename F, typename Diff>
T Newton(F fun, Diff diff, T a, T b, T eps, size_t& iterCount, int method) {
	// Начальная итерация методом хорд
	T fa = fun(a);
	T fb = fun(b);
	T x0 = (fa * b - fb * a) / (fa - fb);
	++iterCount;
	switch (method)
	{
	case 0:		return Newton(fun, diff, a, b, x0, eps, iterCount);//классика
	case 1: 	return NewtonModification(fun, diff, a, b, x0, eps, iterCount);// модификация
	default:
        return {};
	}
}


template<typename T, typename F, typename Diff>
vector<T> Newton(F fun, Diff diff, vector<pair<T, T>> intervals, T eps, int method) {
	vector<T> res;
	res.reserve(intervals.size());
	size_t iterCount = 0;
	for (pair<T, T>& interval : intervals) {
		T a = interval.first;
		T b = interval.second;
		T x = Newton(fun, diff, a, b, eps, iterCount, method);
		res.push_back(x);
	}
	cout << "iterations: " << iterCount << endl;
	return res;
}



// Классический метод Ньютона для системы уравенений
template<typename T, typename Fun, typename J>
vector<T> NewtonSystem(Fun f1, Fun f2, J jacobi, T eps, const vector<T>& x0)
{
	vector<T> xk = x0;
	T coef;

	size_t iterCount = 0;

	do {
	
		matrix<T> jacobiInv = inv(jacobi(xk[0], xk[1]));
		vector<T> F{ f1(xk[0], xk[1]), f2(xk[0], xk[1]) };
		vector<T> jacobiInvMultF = mult(jacobiInv, F);
		vector<T> newXk;
		newXk.reserve(x0.size());
		for (size_t i = 0; i < x0.size(); ++i) {
			newXk.push_back(xk[i] - jacobiInvMultF[i]);
		}

		coef = norm(jacobiInvMultF);

		xk = newXk;

		iterCount++;

	} while (coef > eps);

	cout << "iterations: " << iterCount << endl;

	return xk;
}



//численная матрица
template<typename T, typename F>
auto numJacobi(F FunFirst, F FunSecond, T e) {
	return [&FunFirst, &FunSecond, e](T x1, T x2) {
		T f1Val = FunFirst(x1, x2);
		T f2Val = FunSecond(x1, x2);
		return matrix<T>{
			{ (FunFirst(x1 + e, x2) - f1Val) / e, (FunFirst(x1, x2 + e) - f1Val) / e },
			{ (FunSecond(x1 + e, x2) - f2Val) / e, (FunSecond(x1, x2 + e) - f2Val) / e }
		};
	};
}


//  метод хорд
template<typename T, typename F>
T chord(F f, T a, T b, T eps, size_t& iterCount)
{
	T xk = a;
	T xk1 = b;
	T xk2, coef;
	vector<T> X;
	X.push_back(xk);
	X.push_back(xk1);
	do
	{
		++iterCount;
		xk2 = xk - f(xk) * (xk1 - xk) / (f(xk1) - f(xk));
		X.push_back(xk2);
		coef = abs(xk2 - xk1);
		xk = xk1;
		xk1 = xk2;
	} while (coef > eps);
	degreeX(X, 0.2);
	return xk2;
}


int main()
{
	typedef double Type;

	//vector<Type> x0test4{ -1,-1 };
     vector<Type> x0test4{ 1,1 };
	//vector<Type> x0test4{ -2.433, -9.654 };
	vector<Type> x0test5{ -1., 4. };


	//локализация
	auto intervals = localization(test2<Type>, -1.0, 10.0, 30);
	cout << "Localization:" << endl;
	cout << intervals << endl;
	cout << endl;

	size_t iterCount;
	cout << "Method chord:" << endl;
	auto f = test3<Type>;
	Type rootsChord = chord(f, 0.0, 0.1, 1e-3, iterCount);
	//for (Type& root : rootsChord)
		cout << "x = " << rootsChord << "\n";

	cout << endl;


	cout << "Method bisection:" << endl;
	vector<Type> rootsBis = bisection(test1<Type>, intervals, 1e-6);
	for (Type& root : rootsBis)
		cout << "x = " << root << "\n";

	cout << endl;

	//// первый тест 
	//cout << "Method Newton:" << endl;
 //  // vector<Type> rootsNew = Newton(test2<Type>, numDeriv(test2<Type>, 1e-6), intervals, 1e-6, 0);
 //   vector<Type> rootsNew = Newton(test2<Type>, test2Df<Type>, intervals, 1e-6, 0);
	//for (Type& root : rootsNew)
	//	cout << "x = " << root << "\n";


	//cout << endl;


	// второй тест 
	cout << "Method Newton:" << endl;
	vector<Type> rootsNew = Newton(test2<Type>, numDeriv(test2<Type>, 1e-3), intervals, 1e-6, 0);
	//vector<Type> rootsNew = Newton(test2<Type>, test2Df<Type>, intervals, 1e-6, 0);
	for (Type& root : rootsNew)
		cout << "x = " << root << "\n";


	cout << endl;

    iterCount = 0;
	cout << "x0 = 8, test 2, Newton classic vs modification:" << endl;
	Type rootNewTest2 = Newton(test2<Type>, test2Df<Type>, -1., 10., 8., 1e-6, iterCount);
	cout << "x = " << rootNewTest2 << endl;
    iterCount = 0;
	Type rootNewTest2Mod = NewtonModification(test2<Type>, test2Df<Type>, -1., 10., 8., 1e-6, iterCount);
	cout << "x = " << rootNewTest2Mod << endl;

	cout << endl;

	iterCount = 0;
	cout << "x0 = 0, test 3, Newton classic vs modification:" << endl;
	// тут зацикливание
	//Type rootNewTest3 = Newton(test3<Type>, test3Df<Type>, 0., 1., 0., 1e-3, iterCount);
	//cout << "x = " << rootNewTest3 << endl;	
	// а тут всё ок
	iterCount = 0;
	Type rootNewTest3Mod = NewtonModification(test3<Type>, test3Df<Type>, 0., 1., 0., 1e-3, iterCount);
	cout << "x = " << rootNewTest3Mod << endl;

	cout << endl;

	cout << "Method Newton m = 2:" << endl;
	iterCount = 0;
	Type t0NewAn = Newton(test0<Type>, test0Df<Type>, 0.0, 2.0, 0.2, 1e-6, iterCount, 2);
	cout << "iterations: " << iterCount << endl;
	cout << "x=: " << t0NewAn << endl;
	

	cout << "Method Newton test 4" << endl;
	vector<Type> rootNewSys4 = NewtonSystem(test4f1<Type>, test4f2<Type>,
		//numJacobi<Type>(test4f1<Type>, test4f2<Type>, 1e-6),
		test4Jacobi<Type>,
		1e-6,
		x0test4);
	print_vec(rootNewSys4);

	cout << endl;

	cout << "Method Newton test 5" << endl;
	vector<Type> rootNewSys5 = NewtonSystem(test5f1<Type>, test5f2<Type>,
		numJacobi<Type>(test5f1<Type>, test5f2<Type>, 1e-6),
		//test5Jacobi<Type>,
		1e-6,
		x0test5);
	print_vec(rootNewSys5);


}

