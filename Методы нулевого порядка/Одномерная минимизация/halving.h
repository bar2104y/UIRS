#pragma once


float halving(float a, float b, float e, func f)
{
	double x, y, z, l;
	while (abs(b - a) >= e) 
	{
		x = (a + b) / 2;
		l = abs(b - a);
		y = a + l / 4;
		z = b - l / 4;
		if (f(y) < f(x))
		{
			b = x;
			x = y;
		}
		else {
			if (f(z) < f(x))
			{
				a = x;
				x = z;
			}
			else
			{
				a = y;
				b = z;
			}
		}

	}
	return x;
}