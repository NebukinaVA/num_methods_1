#include "RK3.h"
#include <iostream>

int main()
{
	int i;
	VC equation(0, 1, 2.5, 1.3, 4.2, 2.3, 0.001, 500, 1e-3, 1.0, 1e-5);
	equation.calculate_w_error();
	std::cout << "RK3 METHOD\n" << equation << std::endl;
	std::cin >> i;
	return 0;
}