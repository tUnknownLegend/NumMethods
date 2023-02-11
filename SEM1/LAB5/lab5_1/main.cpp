#include "Source.cpp"
#include "Source1.cpp"

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
		if ((iterCount + 1) % 9 != 0)
			newXk = xk - fun(xk) / diff(xk);
		else {
			// Защита от зацикливания (просто раз в 9 итераций применяем
			// другой метод, в данном случае метод бисекции)
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
		break;
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


int main()
{
	typedef double Type;


	/*vector<Type> x0test4{ -2.433, -9.654 };
	vector<Type> x0test5{ -1., 4. };
	vector<Type> x0test6{ 0.9, 0.5 };*/


	//локализация
	auto intervals = localization(test2<Type>, -1.0, 10.0, 15);
	cout << "Localization:" << endl;
	cout << intervals << endl;
	cout << endl;


	cout << "Method bisection:" << endl;
	vector<Type> rootsBis = bisection(test2<Type>, intervals, 1e-6);
	for (Type& root : rootsBis)
		cout << "x = " << root << "\n";

	cout << endl;

	//// первый тест 
	//cout << "Method Newton:" << endl;
	//vector<Type> rootsNew = Newton(test1<Type>, numDeriv(test1<Type>, 1e-6), intervals, 1e-6, 0);
	////vector<Type> rootsNew = Newton(test1<Type>, test1Df<Type>, intervals, 1e-6, 0);
	//for (Type& root : rootsNew)
	//	cout << "x = " << root << "\n";


	//cout << endl;


	// второй тест 
	cout << "Method Newton:" << endl;
	//vector<Type> rootsNew = Newton(test2<Type>, numDeriv(test2<Type>, 1e-6), intervals, 1e-6, 1);
	vector<Type> rootsNew = Newton(test2<Type>, test2Df<Type>, intervals, 1e-6, 0);
	for (Type& root : rootsNew)
		cout << "x = " << root << "\n";


	cout << endl;

	size_t iterCount = 0;
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
	cout << endl;
	



}

