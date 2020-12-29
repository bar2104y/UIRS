#pragma once
#include <time.h>
#include <random>

/*
float get_random_val(float a, float b)
{
	int s = rand();
	return (float)( (s % abs((int)(b - a)*1000))/1000+a);
}
*/

float sqrt_extr(float a, float b, float e, func f)
{
	float x1, x2, x3, dx, Fmin, Xmin, polinom_min;
	x1 = (b - a) / 3;
	dx = (b - a) / 4;
	x2 = x1 + dx;

	if (f(x1) > f(x2)) x3 = x1 + dx;
	else x3 = x1 - dx;

	do
	{
		if (f(x1) < f(x2) && f(x1) < f(x3))
		{
			Fmin = f(x1); Xmin = x1;
		}
		else if (f(x2) < f(x1) && f(x2) < f(x3))
		{
			Fmin = f(x2); Xmin = x2;
		}
		else
		{
			Fmin = f(x3); Xmin = x3;
		}

		polinom_min = 0.5 * ((x2*x2 - x3*x3) * f(x1) + (x3*x3 - x1*x1) * f(x2) + (x1*x1 - x2*x2) * f(x3)) / ((x2-x3) * f(x1) + (x3-x1) * f(x2) + (x1-x2) * f(x3));
		if (!((Fmin - f(polinom_min)) / f(polinom_min) < e && (Xmin - polinom_min) / polinom_min < e))
			if (polinom_min >= x1 && polinom_min <= x3)
			{
				x2 = polinom_min;
				x1 = x2 - dx / 2;
				x3 = x2 + dx / 2;
			}
			else
			{
				x1 = polinom_min;
				x2 = x1 + dx/2;

				if (f(x1) > f(x2)) x3 = x1 + dx;
				else x3 = x1 - dx;
			}
	} while ((Fmin - f(polinom_min)) / f(polinom_min) < e && (Xmin - polinom_min) / polinom_min < e);

	return polinom_min;
}