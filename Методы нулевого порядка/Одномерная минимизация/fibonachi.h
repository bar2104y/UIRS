#pragma once

typedef float (*func)(float x);

float fib(int i)
{
	if (i < 1) return 0;
	if (i == 1) return 1;

	return fib(i - 1) + fib(i - 2);
}


float fibonachi_extr(float a, float b, float e, func f)
{
	// ОТрезок [a,b]
	// l == e - длина конечного интервала и константа различимости
	// f(x) - функция
	float l = e;
	float l0 = b - a;

	// СЧитать числа Фибоначчи будем прямо в цикле
	int f0 = 1, f1 = 1;
	int n = 1;

	float r, y, z;

	do
	{
		n++;
		r = f0 + f1;
		f0 = f1;
		f1 = r;
	} while (!(r >= l0 / l));

	y = a + fib(n - 2) / fib(n) * (b - a);
	z = a + fib(n - 1) / fib(n) * (b - a);
	for (int k = 0; k < n; k++)
	{
		if (f(y) <= f(z))
		{
			b = z;
			z = y;
			y = a + fib(n - k - 3) / fib(n - k - 1) * (b - a);
		}
		else
		{
			a = y;
			y = z;
			z = a + fib(n - k - 2) / fib(n - k - 1) * (b - a);
		}

		if (k == n - 3)
		{
			float tmp = y;
			y = z;
			z = y;

			if (f(y) <= f(z))b = z;
			else a = y;

			return ((a + b) / 2);
		}
	}
}