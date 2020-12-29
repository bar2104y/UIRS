#pragma once

float dih_extr(float a, float b, float e, func f)
{
	float l = e;
	float y, z;

	while (true)
	{
		y = (a + b - e) / 2;
		z = (a + b + e) / 2;

		if (f(y) <= f(z)) b = z;
		else a = y;

		if (abs(b - a) < 2*l)
			return ((b + a) / 2);
	}

	return ((a + b) / 2);
}