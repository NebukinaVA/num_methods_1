#pragma once
#include <iostream>
#include <cmath>
#include <vector>

// L, R - given constant conditions
// E0, w - variable parameters
// dI/dx = sin(w*x)*E0/L - R*I/L = f(I, x)
// calculation stops when n steps are done

class VC // variable current
{
private:
	double E0, w, L, R; // given conditions (may be variable)
	double I0, x0;
	double h, eps, xmax, prec; 
	int n; // number of steps
	int N = 0; // total steps
	double Xn, Vn;
	int inc = 0;
	int dec = 0; // total inc dec
	int Smin = 0, Smax = 0; // number of string with Smin, Smax
	int hmax, hmin; //number of string with hmax, hmin
	std::vector<double> arg; //x
	std::vector<double> res; //I
	std::vector<double> reshalf; //I with cap
	std::vector<double> steps; //h
	std::vector<double> ss; //S
	std::vector<double> exres; //exact result
	std::vector<int> hinc;  // total step increases
	std::vector<int> hdec;  // total step decreases
	double func(double x, double I)
	{
		return (sin(w*x)*E0 / L - R * I / L);
	}
public:
	VC(double _x0, double _I0, double _L, double _R, double _E0, double _w, double _h, int _n, double _eps, double _xmax, double _prec)
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
	double RK3(double xn, double In, double h)
	{

		double k1 = func(xn, In);
		double k2 = func(xn + h / 2.0, In + h * k1 / 2.0);
		double k3 = func(xn + h, In + h * (2.0 * k2 - k1));
		In += h * (k1 + 4.0 * k2 + k3) / 6.0;
		return In;
	}
	std::vector<double> calculate()
	{
		exres.push_back(I0);
		arg.push_back(x0);
		res.push_back(I0);
		hinc.push_back(0);
		ss.push_back(0.0);
		hdec.push_back(0);
		steps.push_back(h);
		reshalf.push_back(0.0);
		double xn = x0;
		double In = I0;
		int i = 0;
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
				hinc.insert(hinc.begin() + i + 1, inc);
				hdec.insert(hdec.begin() + i + 1, dec);
				reshalf.insert(reshalf.begin() + i + 1, 0.0);
				ss.insert(ss.begin() + i + 1, 0.0);
				i++;
			}
		}
		N = i;
		Xn = arg[i];
		Vn = res[i];
		hmin = i;
		hmax = 1;
		return res;
	}
	std::vector<double> calculate_w_error()
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
		double xn = x0;
		double In = I0;
		double xhalf = x0;
		double Ihalf, reswcap, Inext;
		double S, Sabs;
		int i = 0;
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
				Sabs = abs(S);
				if (i == 0)
				{
					Smin = i + 1;
					Smax = i + 1;
					hmin = i + 1;
					hmax = i + 1;
				}
				else if (Sabs < ss[Smin])
					Smin = i;
				else if (Sabs > ss[Smax])
					Smax = i;
				if ((Sabs >= (eps / 16.0)) && (Sabs <= eps))
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
						hinc.insert(hinc.begin() + i + 2, inc);
						hdec.insert(hdec.begin() + i + 2, ++dec);
					}
					else
					{
						In = RK3(xn, In, h);
					//	xhalf = xn + h / 2.0;
						xn += h;
						hinc.insert(hinc.begin() + i + 2, inc);
						hdec.insert(hdec.begin() + i + 2, dec);
					}

				}
				else if (Sabs < (eps / 16.0))
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
						hinc.insert(hinc.begin() + i + 2, inc);
						hdec.insert(hdec.begin() + i + 2, ++dec);
					}
					else {
						In = RK3(xn, In, h);
					//	xhalf = xn + h / 2.0;
						xn += h;
						h *= 2.0;
						hinc.insert(hinc.begin() + i + 2, ++inc);
						hdec.insert(hdec.begin() + i + 2, dec);
					}
				}
				else if (Sabs > eps)
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
						hinc.insert(hinc.begin() + i + 2, inc);
						hdec.insert(hdec.begin() + i + 2, ++dec);
					}
					else {
						h = h / 2.0;
						In = RK3(xn, In, h);
					//	xhalf = xn + h / 2.0;
						xn += h;
						hinc.insert(hinc.begin() + i + 2, inc);
						hdec.insert(hdec.begin() + i + 2, ++dec);
					}
				}
				i++;
				S *= 8.0;
				if (h < steps[hmin])
					hmin = i + 1;
				if (h > steps[hmax])
					hmax = i + 1;;
				reshalf.insert(reshalf.begin() + i, reswcap);
				ss.insert(ss.begin() + i, S);
				arg.insert(arg.begin() + i, xn);
				res.insert(res.begin() + i, In);
				steps.insert(steps.begin() + i + 1, h);
				exres.insert(exres.begin() + i, ExactSolution(xn));
			}		
		}
		N = i;
		Xn = arg[i];
		Vn = res[i];
		inc = hinc[i];
		dec = hdec[i];
		return res;
	}
	double ExactSolution(double x)
	{
		return (E0*(R*sin(w*x) - L * w*cos(w*x) + L * w*exp(-R * x / L)) / (pow(L, 2.0) * pow(w, 2.0) + pow(R, 2.0)) + I0 * exp(-R * x / L));
	}
	friend std::ostream & operator<<(std::ostream &out, VC &vc)
	{
		if (vc.res.empty())
			out << "There are no calculated results yet.";
		else {
			out << "n  " << " h n-1  " << "    x    " << "      In        " << "        I^        " << "       I      " << "   I - In    " << "       S*      " << "inc " << "dec" << std::endl;
			double globerr;
			for (int i = 0; i < vc.res.size(); i++)
			{
				globerr = abs(vc.exres[i] - vc.res[i]);
				out << i << "  " << vc.steps[i] << "      " << vc.arg[i] << "    " << vc.res[i] << "    " << vc.reshalf[i] << "    " << vc.exres[i]  << "    " << globerr << "   " << vc.ss[i] << "    " << vc.hinc[i] << "  " << vc.hdec[i] << std::endl;
			}
		}
		return out;
	}
};