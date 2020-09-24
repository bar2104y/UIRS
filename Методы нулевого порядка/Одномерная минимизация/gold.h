#pragma once

#define gs 0.381966

float gold_extr(float a, float b, float e, func f)
{
	float y, z;

	y = a + gs * (b - a);
	z = a + b - y;

	do
	{
		if (f(y) <= f(z))
		{
			b = z;
			z = y;
			y = a + b - y;
		}
		else
		{
			a = y;
			y = z;
			z = a + b - z;
		}
	} while (abs(a - b) > e);
	
	return ((a + b) / 2);
}