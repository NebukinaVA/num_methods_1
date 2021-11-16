#include "RK3.h"
#include <iostream>

int main()
{
	/*
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
	*/
	int i;
	VC equation(0, 0, 0.5, 60.0, 4.0, 6.0, 0.001, 5000, 1e-5, 1.0, 1e-4);
	equation.calculate_w_error();
//	equation.calculate();
	std::cout << "RK3 METHOD\n" << equation << std::endl;
	std::cin >> i;
	return 0;
}