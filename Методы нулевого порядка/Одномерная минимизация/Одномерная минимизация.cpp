#include <iostream>
#include "fibonachi.h"


using namespace std;

//Уравнение параболы, смещенной влево на 1 и вверх на 1
float parabola(float x)
{
	return (x * x + 2 * x + 2); // (x+1)^2 + 1
}

int main()
{
	float e = fibonachi_extr(-10, 10, 0.001, *parabola);
	cout << "x=" << e << "		y=" << parabola(e);
	return 0;
}