#pragma once


/// <summary>
/// Метод равномерного поиска
/// </summary>
/// <param name="a">Начало отрезка</param>
/// <param name="b">Конец отрезка</param>
/// <param name="e">Точность/шаг поиска</param>
/// <param name="f">Ссылка на функцию одного аргумента</param>
/// <returns>Xmin</returns>
float evenly(float a, float b, float e, func f)
{

	float x = a, fx = f(a),	// Текущие значения
		xmin = a, fmin = f(a); // Минимальные значения
	while (x < b)
	{
		x += e;
		fx = f(x);

		if (fx < fmin)
		{
			xmin = x;
			fmin = fx;
		}
	}
	return xmin;
}