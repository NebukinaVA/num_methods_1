#pragma once
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <iterator>

// L, R - given constant conditions
// E0, w - variable parameters
// dI/dx = sin(w*x)*E0/L - R*I/L = f(I, x)
// calculation stops when n steps are done

class VC // variable current
{
private:
	long double E0, w, L, R; // given conditions (may be variable)
	long double I0, x0;
	long double h, eps, xmax; 
	long int n; // number of steps
	std::vector<long double> arg; //x
	std::vector<long double> res; //I
	std::vector<long double> steps; //h
	std::vector<long double> ss; //S
	std::vector<long double> exres; //exact result
	long double func(long double x, long double I)
	{
		return (sin(w*x)*E0 / L - R * I / L);
	}
public:
	VC(long double _x0, long double _I0, long double _L, long double _R, long double _E0, long double _w, long double _h, long int _n, long double _eps, long double _xmax)
	{
		x0 = _x0;
		I0 = _I0;
		L = _L;
		R = _R;
		E0 = _E0;
		w = _w;
		h = _h;
		n = _n;
		eps = _eps;
		xmax = _xmax;
	}
	long double RK3(long double xn, long double In, long double h)
	{

		long double k1 = func(xn, In);
		long double k2 = func(xn + h / 2.0, In + h * k1 / 2.0);
		long double k3 = func(xn + h, In + h * (-k1 + 2.0 * k2));
		In += h * (k1 + 4.0 * k2 + k3) / 6.0;
		return In;
	}
	std::vector<long double> calculate()
	{
		exres.push_back(I0);
		steps.push_back(h);
		arg.push_back(x0);
		res.push_back(I0);
		long double exI = I0;
		long double xn = x0;
		long double In = I0;
		for (long int i = 0; i < n; i++)
		{
			In = RK3(xn, In, h);
			xn += h;
			arg.insert(arg.begin() + i + 1, xn);
			res.insert(res.begin() + i + 1, In);
			steps.insert(steps.begin() + i + 1, h);
			exI = ExactSolution(xn, exI);
			exres.insert(exres.begin() + i + 1, exI);
		}
		return res;
	}
	std::vector<long double> calculate_w_error()
	{
		steps.push_back(0.0);
		arg.push_back(x0);
		res.push_back(I0);
		long double exI = I0;
		long double xn = x0;
		long double In = I0;
		long double xhalf = x0;
		long double Ihalf;
		long double S;
		for (long int i = 0; i < n; i++)
		{
			Ihalf = RK3(xn, In, h / 2.0);
			S = (RK3(xhalf, Ihalf, h / 2.0) - RK3(xn, In, h)) / 7.0;
			steps.insert(steps.begin() + i + 1, S);
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
			arg.insert(arg.begin() + i + 1, xn);
			res.insert(res.begin() + i + 1, In);
			steps.insert(steps.begin() + i + 1, h);
			exI = ExactSolution(xn, exI);
			exres.insert(exres.begin() + i + 1, exI);
		}
		return res;
	}
	long double ExactSolution(long double x0, long double I0)
	{
		return ((E0*exp(R / L)*(sin(w*x0)*R - w * L*cos(w*x0)) + E0 * w*L) / (pow(w, 2.0) * L * 2 + pow(R, 2.0)) + I0);
	}
	friend std::ostream & operator<<(std::ostream &out, VC &vc)
	{
		if (vc.res.empty())
			out << "There are no calculated results yet.";
		else {
			for (long int i = 0; i < vc.res.size(); i++)
			{
				out << "x" << i << " = " << vc.arg[i] << ",  " << "I" << i << " = " << vc.res[i] << std::endl;
			}
		}
		return out;
	}
};