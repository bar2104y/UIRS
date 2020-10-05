#pragma once



// Fn - n-ое число фибоначи
float fib(int i)
{
	if (i < 1) return 0;	// F0 = 0
	if (i == 1) return 1;	// F0 = 1

	// F(n) = F(n-1) + F(n-2)
	return fib(i - 1) + fib(i - 2);
}


float fibonachi_extr(float a, float b, float e, func f)
{
	// ОТрезок [a,b]
	// l == e - длина конечного интервала и константа различимости
	// f(x) - функция
	float l = e;
	float l0 = b - a;

	// Считать числа Фибоначчи будем прямо в цикле
	int f0 = 1, f1 = 1;
	int n = 1;


	float r, y, z;

	// Найдем n - число повторений вычислений
	do
	{
		n++;
		r = f0 + f1;
		f0 = f1;
		f1 = r;
	} while (!(r >= l0 / l));


	y = a + fib(n - 2) / fib(n) * (b - a); // Новая левая граница
	z = a + fib(n - 1) / fib(n) * (b - a); // Новая правая граница
	
	
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

		// Условие выхода
		if (k == n - 3)
		{
			float tmp = y;
			y = z;
			z = tmp+e;

			
			if (f(y) <= f(z)) b = z;
			else a = y;

			return ((a + b) / 2);
		}
	}
}