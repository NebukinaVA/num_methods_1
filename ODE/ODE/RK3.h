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
	long double h, eps, xmax, prec; 
	long int n; // number of steps
	std::vector<long double> arg; //x
	std::vector<long double> res; //I
	std::vector<long double> reshalf; //I with cap
	std::vector<long double> steps; //h
	std::vector<long double> ss; //S
	std::vector<long double> exres; //exact result
	std::vector<long int> hinc;  // total step increases
	std::vector<long int> hdec;  // total step decreases
	long double func(long double x, long double I)
	{
		return (sin(w*x)*E0 / L - R * I / L);
	}
public:
	VC(long double _x0, long double _I0, long double _L, long double _R, long double _E0, long double _w, long double _h, long int _n, long double _eps, long double _xmax, long double _prec)
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
		prec = _prec;
	}
	long double RK3(long double xn, long double In, long double h)
	{

		long double k1 = func(xn, In);
		long double k2 = func(xn + h / 2.0, In + h * k1 / 2.0);
//		long double k3 = func(xn + h, In + h * (-k1 + 2.0 * k2));
		long double k3 = func(xn + h, In + h * (2.0 * k2 - k1));
		In += h * (k1 + 4.0 * k2 + k3) / 6.0;
		return In;
	}
	std::vector<long double> calculate()
	{
		exres.push_back(I0);
		arg.push_back(x0);
		res.push_back(I0);
		hinc.push_back(0);
		ss.push_back(0.0);
		hdec.push_back(0);
		steps.push_back(h);
		reshalf.push_back(0.0);
		long double xn = x0;
		long double In = I0;
		long int i = 0;
		while (i < n)
		{
			if ((xn > (xmax - prec)) && (xn < xmax))
			{
				break;
			}
			else {
				if ((xn + h) > xmax)
				{
					while (((xn + h) > xmax) && (xn < (xmax - prec)))
					{
						h /= 2.0;
					}
					In = RK3(xn, In, h);
					xn += h;
				}
				else
				{
					In = RK3(xn, In, h);
					xn += h;
				}
				arg.insert(arg.begin() + i + 1, xn);
				res.insert(res.begin() + i + 1, In);
				steps.insert(steps.begin() + i + 1, h);
				exres.insert(exres.begin() + i + 1, ExactSolution(xn));
				hinc.insert(hinc.begin() + i + 1, 0);
				hdec.insert(hdec.begin() + i + 1, 0);
				reshalf.insert(reshalf.begin() + i + 1, 0.0);
				ss.insert(ss.begin() + i + 1, 0.0);
				i++;
			}
		}
		return res;
	}
	std::vector<long double> calculate_w_error()
	{
		exres.push_back(I0);
		ss.push_back(0.0);
		arg.push_back(x0);
		res.push_back(I0);
		steps.push_back(0.0);
		steps.push_back(h);
		reshalf.push_back(0.0);
		hinc.push_back(0);
		hdec.push_back(0);
		hinc.push_back(0);
		hdec.push_back(0);
		long double xn = x0;
		long double In = I0;
		long double xhalf = x0;
		long double Ihalf, reswcap, Inext;
		long double S;
		long int i = 0;
		while (i < n)
		{
			if ((xn > (xmax - prec)) && (xn < xmax))
			{
				break;
			}			
			else
			{
				Ihalf = RK3(xn, In, h / 2.0);
				xhalf = xn + h / 2.0;
				reswcap = RK3(xhalf, Ihalf, h / 2.0);
				Inext = RK3(xn, In, h);
				S = (reswcap - Inext) / 7.0;
				if ((abs(S) >= (eps / 16.0)) && (abs(S) <= eps))
				{
					if ((xn + h) > xmax)
					{
						while (((xn + h) > xmax) && (xn < (xmax - prec)))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						In = RK3(xn, In, h);
						hinc.insert(hinc.begin() + i + 2, hinc[i+1]);
						hdec.insert(hdec.begin() + i + 2, ++(hdec[i+1]));
					}
					else
					{
						In = RK3(xn, In, h);
					//	xhalf = xn + h / 2.0;
						xn += h;
						hinc.insert(hinc.begin() + i + 2, hinc[i+1]);
						hdec.insert(hdec.begin() + i + 2, hdec[i+1]);
					}

				}
				else if (abs(S) < (eps / 16.0))
				{
					if ((xn + h) > xmax)
					{
						while (((xn + h) > xmax) && (xn < (xmax - prec)))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						In = RK3(xn, In, h);
						hinc.insert(hinc.begin() + i + 2, hinc[i+1]);
						hdec.insert(hdec.begin() + i + 2, ++(hdec[i+1]));
					}
					else {
						In = RK3(xn, In, h);
					//	xhalf = xn + h / 2.0;
						xn += h;
						h *= 2.0;
						hinc.insert(hinc.begin() + i + 2, ++(hinc[i+1]));
						hdec.insert(hdec.begin() + i + 2, hdec[i+1]);
					}
				}
				else if (abs(S) > eps)
				{
					if ((xn + h) > xmax)
					{
						while (((xn + h) > xmax) && (xn < (xmax - prec)))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						In = RK3(xn, In, h);
						hinc.insert(hinc.begin() + i + 2, hinc[i+1]);
						hdec.insert(hdec.begin() + i + 2, ++(hdec[i+1]));
					}
					else {
						h = h / 2.0;
						In = RK3(xn, In, h);
					//	xhalf = xn + h / 2.0;
						xn += h;
						hinc.insert(hinc.begin() + i + 2, hinc[i+1]);
						hdec.insert(hdec.begin() + i + 2, ++(hdec[i+1]));
					}
				}
				i++;
				S *= 8.0;
				reshalf.insert(reshalf.begin() + i, reswcap);
				ss.insert(ss.begin() + i, S);
				arg.insert(arg.begin() + i, xn);
				res.insert(res.begin() + i, In);
				steps.insert(steps.begin() + i + 1, h);
				exres.insert(exres.begin() + i, ExactSolution(xn));
			}		
		}
		return res;
	}
	long double ExactSolution(long double x)
	{
		return (E0*(R*sin(w*x) - L * w*cos(w*x) + L * w*exp(-R * x / L)) / (pow(L, 2.0) * pow(w, 2.0) + pow(R, 2.0)) + I0 * exp(-R * x / L));
	}
	friend std::ostream & operator<<(std::ostream &out, VC &vc)
	{
		if (vc.res.empty())
			out << "There are no calculated results yet.";
		else {
			out << "n  " << " h n-1  " << "    x    " << "      In        " << "        I^        " << "       I      " << "   I - In    " << "       S*      " << "inc " << "dec" << std::endl;
			long double globerr;
			for (long int i = 0; i < vc.res.size(); i++)
			{
				globerr = abs(vc.exres[i] - vc.res[i]);
				out << i << "  " << vc.steps[i] << "      " << vc.arg[i] << "    " << vc.res[i] << "    " << vc.reshalf[i] << "    " << vc.exres[i]  << "    " << globerr << "   " << vc.ss[i] << "    " << vc.hinc[i] << "  " << vc.hdec[i] << std::endl;
			}
		}
		return out;
	}
};