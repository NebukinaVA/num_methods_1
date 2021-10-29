#include "RK3.h"
#include <iostream>

int main()
{
	int i;
	VC equation(500);
	equation.calculate(0.01);
	std::cout << "RK3 METHOD\n" << equation << std::endl;
	std::cin >> i;
	return 0;
}