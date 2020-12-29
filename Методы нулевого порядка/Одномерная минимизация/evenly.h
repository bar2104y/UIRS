#pragma once


/// <summary>
/// ����� ������������ ������
/// </summary>
/// <param name="a">������ �������</param>
/// <param name="b">����� �������</param>
/// <param name="e">��������/��� ������</param>
/// <param name="f">������ �� ������� ������ ���������</param>
/// <returns>Xmin</returns>
float evenly(float a, float b, float e, func f)
{

	float x = a, fx = f(a),	// ������� ��������
		xmin = a, fmin = f(a); // ����������� ��������
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