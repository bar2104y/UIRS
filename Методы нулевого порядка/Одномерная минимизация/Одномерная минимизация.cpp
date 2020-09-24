#include <iostream>
using namespace std;


typedef float (*func)(float x);

#include "fibonachi.h";
#include "dihotomia.h";
#include "gold.h";




//Уравнение параболы, смещенной влево на 1 и вверх на 1
float parabola(float x)
{
	return (x * x + 2 * x + 2); // (x+1)^2 + 1
}

int main()
{
	float e;
	e = fibonachi_extr(-10, 10, 0.001, *parabola);
	cout << "x=" << e << "		y=" << parabola(e) << endl;

	e = dih_extr(-10, 10, 0.001, *parabola);
	cout << "x=" << e << "		y=" << parabola(e) << endl;

	e = gold_extr(-10, 10, 0.001, *parabola);
	cout << "x=" << e << "		y=" << parabola(e) << endl;
	
	return 0;
}