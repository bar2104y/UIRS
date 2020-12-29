#include <iostream>
using namespace std;


typedef float (*func)(float x);


#include "fibonachi.h"
#include "dihotomia.h"
#include "gold.h"
#include "kv_intr.h"
#include "evenly.h"
#include "halving.h"


/// +Равномерный
/// +Делением полоам
/// +Дихотомия
/// +Золотое сечение
/// +Фибоначи



//Присер из пособия
float parabola(float x)
{
	return (2*x * x - 12 * x);
}



int main()
{
	float e;
	e = fibonachi_extr(-5, 10, 0.001, *parabola);
	cout << "FIBONACHI	:	x=" << e << "		y=" << parabola(e)<< endl;

	e = dih_extr(0, 10, 0.001, *parabola);
	cout << "DIHOTOMIYA	:	x=" << e << "		y=" << parabola(e) << endl;

	e = gold_extr(-5, 10, 0.001, *parabola);
	cout << "ZOL SECH	:	x=" << e << "		y=" << parabola(e) << endl;

	e = evenly(-5, 10, 0.001, *parabola);
	cout << "RAVNOMERNY	:	x=" << e << "		y=" << parabola(e) << endl;

	e = halving(-5, 10, 0.001, *parabola);
	cout << "DEL POPOLAM	:	x=" << e << "		y=" << parabola(e) << endl;

	e = sqrt_extr(-5, 10, 0.001, *parabola);
	cout << "KVADR INTER	:	x=" << e << "		y=" << parabola(e) << endl;
	
	return 0;
}