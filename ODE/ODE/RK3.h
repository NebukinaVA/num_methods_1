#pragma once
#include <iostream>
#include <cmath>
#include <map>
#include <iterator>

// L, R - по условию заданы и постоянны 
// E0, w - изменяемые параметры
// dI/dx = sin(w*x)*E0/L - R*I/L = f(I, x)
// остановка счета по выполнению заданного числа шагов

class VC // variable current
{
private:
	double E0, w;
	double L = 2.5;
	double R = 6.21;
	long int n; // кол-во шагов
	std::map<long double, long double> results;  // (x, I)
	long double func(long double x, long double I)
	{
		return (sin(w*x)*E0 / L - R * I / L);
	}
public:
	VC(long int _n, double _E0 = 4.2, double _w = 3.2)
	{
		E0 = _E0;
		w = _w;
		n = _n;
	}
	long double RK3(long double h, long double xn = 0.0, long double In = 1.0)
	{

		long double k1 = func(xn, In);
		long double k2 = func(xn + h / 2.0, In + h * k1 / 2.0);
		long double k3 = func(xn + h, In + h * (-k1 + 2.0 * k2));
		In += h * (k1 + 4.0 * k2 + k3) / 6.0;
		return In;
	}
	std::map<long double, long double> calculate(long double h, long double eps = 1e-3, long double x0 = 0.0, long double I0 = 1.0)
	{
		results.insert(std::make_pair(x0, I0));
		long double xn = x0;
		long double In = I0;
		long double xhalf = x0;
		long double Ihalf;
		long double S;
		for (int i = 0; i < n; i++)
		{
			Ihalf = RK3(xn, In, h / 2.0);
			S = (RK3(xhalf, Ihalf, h / 2.0) - RK3(xn, In, h)) / 7.0;
			if ((abs(S) >= (eps / 16.0)) && (abs(S) <= eps))
			{
				In = RK3(xn, In, h);
			}
			else if (abs(S) < (eps / 16.0))
			{
				In = RK3(xn, In, h);
				h *= 2.0;
			}
			else if (abs(S) > eps)
			{
				h = h / 2.0;
				In = RK3(xn, In, h);
			}
			xhalf = xn + h / 0.2;
			xn += h;
			results.insert(std::make_pair(xn, In));
			
		}
		return results;
	}
	friend std::ostream & operator<<(std::ostream &out, VC &vc)
	{
		if (vc.results.empty())
			out << "There are no calculated results yet.";
		else {
			std::map<long double, long double>::iterator it;
			long int i = 0;
			for (it = vc.results.begin(); it != vc.results.end(); ++it)
			{
				out << "x" << i << " = " << it->first << ",  " << "I" << i << " = " << it->second << std::endl;
				++i;
			}
		}
		return out;
	}
};