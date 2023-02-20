
// Главный файл lab_one

#include "main.h"
#include "diffure.h"
#include "matrix_operations.h"

int main() {

	
	vector<double> u0 = { 1.0, 0.0 };
	int n = 5;
	double t0 = 0.0;
	double T = 5.0;
	double h = (T - t0) / n;
	int dim = 2;
	matrix<double> y(n, vector<double>(dim));
	matrix<double> z(n, vector<double>(dim));
	matrix<double> t(n, vector<double>(dim));
	matrix<double> p(n, vector<double>(dim));
	matrix<double> q(n, vector<double>(dim));
	matrix<double> r(n, vector<double>(dim));
	matrix<double> s(n, vector<double>(dim));

  

	cout << "\n Exp Euler method: \n";
	y = Euler_explicit(u0, n, h, dim);
	print(y);
    //НОВОЕ ВРОДЕ КАК СДЕЛАЛА РАЗНОСТЬ 
    cout << "f-f1 with h: " << endl;
    for (int i = 0; i < y.size(); i++)
    {
       cout  << cos(i * h) - y[i][0] << endl;
   
    }
    
    matrix<double> y2(2*n, vector<double>(dim));
    y2 = Euler_explicit(u0, 2*n, h/2, dim);
    /*for (int i = 0; i < y2.size(); i+=2) {
        cout << i*h/2 << endl;
    }*/
    cout << "f-f2 with h/2: " << endl;
    for (int i = 0; i < y2.size(); i+=2)
    {
       cout << cos(i * h / 2) - y2[i][0] << endl;

    }
    cout << "p: " << endl;
    cout << log2((-y[3][0] + cos(3 * h)) / (-y2[6][0] + cos(6 * h/2))) << endl;

    matrix<double> y3(4 * n, vector<double>(dim));
    y3 = Euler_explicit(u0, 4 * n, h / 4, dim);
    cout << "f-f3 with h/4: " << endl;
    for (int i = 0; i < y3.size(); i += 4)
    {
        cout << cos(i * h / 4) - y3[i][0] << endl;

    }

    cout << "p: " << endl;
    cout << log2((-y2[6][0] + cos(6 * h / 2)) / (-y3[12][0] + cos(12 * h / 4))) << endl;
   
	std::string filename1 = "C:\\Users\\kate\\Documents\\GitHub\\NumMethods\\SEM2\\LABA1\\laba1\\Euler_explicit.txt";
	std::ofstream Euler1;
	Euler1.open(filename1);
	for (int i = 0; i < y.size(); i++)
	{
		Euler1 << "{";
		for (int j = 0; j < y[0].size(); j++)
		{
			Euler1 << y[i][j];
			if (j != y[0].size() - 1) Euler1 << ",";
		}
		Euler1 << "}, ";
	}
	Euler1.close();


	cout << "\n  Imp Euler method:" << endl;
	q = Euler_implicit(u0, n, h, dim);
	print(q);
    //НОВОЕ ВРОДЕ КАК СДЕЛАЛА РАЗНОСТЬ 
    cout << "f-f1 with h: " << endl;
    for (int i = 0; i < q.size(); i++)
    {
        cout << cos(i * h) - q[i][0] << endl;

    }

    matrix<double> q2(2 * n, vector<double>(dim));
    q2 = Euler_implicit(u0, 2 * n, h / 2, dim);
    /*for (int i = 0; i < y2.size(); i+=2) {
        cout << i*h/2 << endl;
    }*/
    cout << "f-f2 with h/2: " << endl;
    for (int i = 0; i < q2.size(); i += 2)
    {
        cout << cos(i * h / 2) - q2[i][0] << endl;

    }
    cout << "p: " << endl;
    cout << log2((-q[3][0] + cos(3 * h)) / (-q2[6][0] + cos(6 * h / 2))) << endl;

    matrix<double> q3(4 * n, vector<double>(dim));
    q3 = Euler_implicit(u0, 4 * n, h / 4, dim);
    cout << "f-f3 with h/4: " << endl;
    for (int i = 0; i < q3.size(); i += 4)
    {
        cout << cos(i * h / 4) - q3[i][0] << endl;

    }

    cout << "p: " << endl;
    cout << log2((-q2[6][0] + cos(6 * h / 2)) / (-q3[12][0] + cos(12 * h / 4))) << endl;
	std::string filename2 = "C:\\Users\\kate\\Documents\\GitHub\\NumMethods\\SEM2\\LABA1\\laba1\\Euler_implicit.txt";
	std::ofstream Euler2;
	Euler2.open(filename2);
	for (int i = 0; i < q.size(); i++)
	{
		Euler2 << "{";
		for (int j = 0; j < q[0].size(); j++)
		{
			Euler2 << q[i][j];
			if (j != q[0].size() - 1) Euler2 << ",";
		}
		Euler2 << "}, ";
	}
	Euler2.close();

	cout << "\n Symmetric:" << endl;
	r = symmetric(u0, h, n, dim);
	print(r);
    //НОВОЕ ВРОДЕ КАК СДЕЛАЛА РАЗНОСТЬ 
    cout << "f-f1 with h: " << endl;
    for (int i = 0; i < r.size(); i++)
    {
        cout << cos(i * h) - r[i][0] << endl;

    }

    matrix<double> r2(2 * n, vector<double>(dim));
    r2 = symmetric(u0, h / 2, 2 * n, dim);
    /*for (int i = 0; i < y2.size(); i+=2) {
        cout << i*h/2 << endl;
    }*/
    cout << "f-f2 with h/2: " << endl;
    for (int i = 0; i < r2.size(); i += 2)
    {
        cout << cos(i * h / 2) - r2[i][0] << endl;

    }
    cout << "p: " << endl;
    cout << log2((-r[3][0] + cos(3 * h)) / (-r2[6][0] + cos(6 * h / 2))) << endl;

    matrix<double> r3(4 * n, vector<double>(dim));
    r3 = symmetric(u0, h / 4, 4 * n, dim);
    cout << "f-f3 with h/4: " << endl;
    for (int i = 0; i < r3.size(); i += 4)
    {
        cout << cos(i * h / 4) - r3[i][0] << endl;

    }

    cout << "p: " << endl;
    cout << log2((-r2[6][0] + cos(6 * h / 2)) / (-r3[12][0] + cos(12 * h / 4))) << endl;
	std::string filename3 = "C:\\Users\\kate\\Documents\\GitHub\\NumMethods\\SEM2\\LABA1\\laba1\\Symmetric.txt";
	std::ofstream sym;
	sym.open(filename3);
	for (int i = 0; i < r.size(); i++)
	{
		sym << "{";
		for (int j = 0; j < r[0].size(); j++)
		{
			sym << r[i][j];
			if (j != r[0].size() - 1) sym << ",";
		}
		sym << "}, ";
	}
	sym.close();

	cout << "\n Method R-K 2: ";
	s = Runge_Kutta2(u0, h, n, dim);
	print(s);
    //НОВОЕ ВРОДЕ КАК СДЕЛАЛА РАЗНОСТЬ 
    cout << "f-f1 with h: " << endl;
    for (int i = 0; i < s.size(); i++)
    {
        cout << cos(i * h) - s[i][0] << endl;

    }

    matrix<double> s2(2 * n, vector<double>(dim));
    s2 = Runge_Kutta2(u0, h / 2, 2 * n, dim);
    /*for (int i = 0; i < y2.size(); i+=2) {
        cout << i*h/2 << endl;
    }*/
    cout << "f-f2 with h/2: " << endl;
    for (int i = 0; i < s2.size(); i += 2)
    {
        cout << cos(i * h / 2) - s2[i][0] << endl;

    }
    cout << "p: " << endl;
    cout << log2((-s[3][0] + cos(3 * h)) / (-s2[6][0] + cos(6 * h / 2))) << endl;

    matrix<double> s3(4 * n, vector<double>(dim));
    s3 = Runge_Kutta2(u0, h / 4, 4 * n, dim);
    cout << "f-f3 with h/4: " << endl;
    for (int i = 0; i < s3.size(); i += 4)
    {
        cout << cos(i * h / 4) - s3[i][0] << endl;

    }

    cout << "p: " << endl;
    cout << log2((-s2[6][0] + cos(6 * h / 2)) / (-s3[12][0] + cos(12 * h / 4))) << endl;
	std::string filename4 = "C:\\Users\\kate\\Documents\\GitHub\\NumMethods\\SEM2\\LABA1\\laba1\\Runge_Kutta2.txt";
	std::ofstream RK2;
	RK2.open(filename4);
	for (int i = 0; i < s.size(); i++)
	{
		RK2 << "{";
		for (int j = 0; j < s[0].size(); j++)
		{
			RK2 << s[i][j];
			if (j != s[0].size() - 1) RK2 << ",";
		}
		RK2 << "}, ";
	}
	RK2.close();

	cout << "\n Method R-K 4: ";
	z = Runge_Kutta4(u0, h, n, dim);
	print(z);
    //НОВОЕ ВРОДЕ КАК СДЕЛАЛА РАЗНОСТЬ 
    cout << "f-f1 with h: " << endl;
    for (int i = 0; i < z.size(); i++)
    {
        cout << cos(i * h) - z[i][0] << endl;

    }

    matrix<double> z2(2 * n, vector<double>(dim));
    z2 = Runge_Kutta4(u0, h / 2, 2 * n, dim);
    /*for (int i = 0; i < y2.size(); i+=2) {
        cout << i*h/2 << endl;
    }*/
    cout << "f-f2 with h/2: " << endl;
    for (int i = 0; i < z2.size(); i += 2)
    {
        cout << cos(i * h / 2) - z2[i][0] << endl;

    }
    cout << "p: " << endl;
    cout << log2((-z[3][0] + cos(3 * h)) / (-z2[6][0] + cos(6 * h / 2))) << endl;

    matrix<double> z3(4 * n, vector<double>(dim));
    z3 = Runge_Kutta4(u0, h / 4, 4 * n, dim);
    cout << "f-f3 with h/4: " << endl;
    for (int i = 0; i < z3.size(); i += 4)
    {
        cout << cos(i * h / 4) - z3[i][0] << endl;

    }

    cout << "p: " << endl;
    cout << log2((-z2[6][0] + cos(6 * h / 2)) / (-z3[12][0] + cos(12 * h / 4))) << endl;
	std::string filename5 = "C:\\Users\\kate\\Documents\\GitHub\\NumMethods\\SEM2\\LABA1\\laba1\\Runge_Kutta4.txt";
	std::ofstream RK4;
	RK4.open(filename5);
	for (int i = 0; i < z.size(); i++)
	{
		RK4 << "{";
		for (int j = 0; j < z[0].size(); j++)
		{
			RK4 << z[i][j];
			if (j != z[0].size() - 1) RK4 << ",";
		}
		RK4 << "}, ";
	}
	RK4.close();

	cout << "\n Method Adams 4:" << endl;
	t = Adams(u0, h, n, dim);
    print(t);
    //НОВОЕ ВРОДЕ КАК СДЕЛАЛА РАЗНОСТЬ 
    cout << "f-f1 with h: " << endl;
    for (int i = 0; i < t.size(); i++)
    {
        cout << cos(i * h) - t[i][0] << endl;

    }

    matrix<double> t2(2 * n, vector<double>(dim));
    t2 = Adams(u0, h / 2, 2 * n, dim);
    cout << "f-f2 with h/2: " << endl;
    for (int i = 0; i < t2.size(); i += 2)
    {
        cout << cos(i * h / 2) - t2[i][0] << endl;

    }
    cout << "p: " << endl;
    cout << log2((-t[3][0] + cos(3 * h)) / (-t2[6][0] + cos(6 * h / 2))) << endl;

    matrix<double> t3(4 * n, vector<double>(dim));
    t3 = Adams(u0, h / 4, 4 * n, dim);
    cout << "f-f3 with h/4: " << endl;
    for (int i = 0; i < t3.size(); i += 4)
    {
        cout << cos(i * h / 4) - t3[i][0] << endl;

    }

    cout << "p: " << endl;
    cout << log2((-t2[6][0] + cos(6 * h / 2)) / (-t3[12][0] + cos(12 * h / 4))) << endl;
	std::string filename6 = "C:\\Users\\kate\\Documents\\GitHub\\NumMethods\\SEM2\\LABA1\\laba1\\Adams.txt";
	std::ofstream adam;
	adam.open(filename6);
	for (int i = 0; i < t.size(); i++)
	{
		adam << "{";
		for (int j = 0; j < t[0].size(); j++)
		{
			adam << t[i][j];
			if (j != t[0].size() - 1) adam << ",";
		}
		adam << "}, ";
	}
	adam.close();

	cout << "\n Method prediction:" << endl;
	p = prediction(u0, h, n, dim);
	print(p);
    //НОВОЕ ВРОДЕ КАК СДЕЛАЛА РАЗНОСТЬ 
    cout << "f-f1 with h: " << endl;
    for (int i = 0; i < p.size(); i++)
    {
        cout << cos(i * h) - p[i][0] << endl;

    }

    matrix<double> p2(2 * n, vector<double>(dim));
    p2 = prediction(u0, h / 2, 2 * n, dim);
    cout << "f-f2 with h/2: " << endl;
    for (int i = 0; i < p2.size(); i += 2)
    {
        cout << cos(i * h / 2) - p2[i][0] << endl;

    }
    cout << "p: " << endl;
    cout << log2((-p[3][0] + cos(3 * h)) / (-p2[6][0] + cos(6 * h / 2))) << endl;

    matrix<double> p3(4 * n, vector<double>(dim));
    p3 = prediction(u0, h / 4, 4 * n, dim);
    cout << "f-f3 with h/4: " << endl;
    for (int i = 0; i < p3.size(); i += 4)
    {
        cout << cos(i * h / 4) - p3[i][0] << endl;

    }

    cout << "p: " << endl;
    cout << log2((-p2[6][0] + cos(6 * h / 2)) / (-p3[12][0] + cos(12 * h / 4))) << endl;
	std::string filename7 = "C:\\Users\\kate\\Documents\\GitHub\\NumMethods\\SEM2\\LABA1\\laba1\\Prediction.txt";
	std::ofstream pred;
	pred.open(filename7);
	for (int i = 0; i < p.size(); i++)
	{
		pred << "{";
		for (int j = 0; j < p[0].size(); j++)
		{
			pred << p[i][j];
			if (j != p[0].size() - 1) pred << ",";
		}
		pred << "}, ";
	}
	pred.close();
}