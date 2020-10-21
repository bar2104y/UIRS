#pragma once
#include <iostream>
#include <iomanip>

using namespace std;

class FloatMatrix
{
public:
	FloatMatrix();
	FloatMatrix(const FloatMatrix &obj);
	FloatMatrix(unsigned int a);
	FloatMatrix(unsigned int a, unsigned int b);
	FloatMatrix(unsigned int a, unsigned int b, float** ref);
	~FloatMatrix();

	void FillFloat(float v);

	void Print();
	void PrintRow(unsigned int a);
	//void PrintColumn(unsigned int n);

	FloatMatrix operator + (FloatMatrix &m2)
	{
		unsigned int maxM, maxN, minM,minN;
		if (this->m > m2.m) {
			maxM = this->m;
			minM = m2.m;
		}
		else
		{
			maxM = m2.m;
			minM = this->m;
		}

		if (this->n > m2.n) {
			maxN = this->n;
			minN = m2.n;
		}
		else
		{
			maxN = m2.n;
			minN = this->n;
		}
		// Создать новую матрицу максимальной размерности и в нее класть результат
		
		FloatMatrix res = FloatMatrix(maxM, maxN);

		for (unsigned int i = 0; i < this->m; i++)
		{
			for (unsigned int j = 0; j < this->n; j++)
				res.matrix[i][j] = this->matrix[i][j];
		}

		for (unsigned int i = 0; i < m2.m; i++)
		{
			for (unsigned int j = 0; j < m2.n; j++)
				res.matrix[i][j] += m2.matrix[i][j];
		}

		
		return FloatMatrix(res.m,res.n, res.matrix);

	}



private:
	unsigned int m, n;//строка, столбец
	float** matrix;
};

FloatMatrix::FloatMatrix()
{
	n = 1;
	m = 1;
	matrix = new float*[m];
	matrix[0] = new float[n];
	matrix[0][0] = 0;
}

FloatMatrix::FloatMatrix(const FloatMatrix &obj)
{
	m = obj.m;
	n = obj.n;
	matrix = new float* [m];
	for (unsigned int i = 0; i < m; i++)
	{
		matrix[i] = new float[n];
		for (unsigned int j = 0; j < n; j++)
			matrix[i][j] = 0;
	}

}

FloatMatrix::FloatMatrix(unsigned int a)
{
	n = a;
	m = a;
	matrix = new float* [m];
	for (unsigned int i = 0; i < m; i++)
	{
		matrix[i] = new float[n];
		for (unsigned int j = 0; j < n; j++)
			matrix[i][j] = 0;
	}
}

FloatMatrix::FloatMatrix(unsigned int a, unsigned int b)
{
	m = a;
	n = b;
	matrix = new float* [m];
	for (unsigned int i = 0; i < m; i++)
	{
		matrix[i] = new float[n];
		for (unsigned int j = 0; j < n; j++)
			matrix[i][j] = 0;
	}
}

FloatMatrix::FloatMatrix(unsigned int a, unsigned int b, float** ref)
{
	m = a;
	n = b;
	matrix = new float* [m];
	for (unsigned int i = 0; i < m; i++)
	{
		matrix[i] = new float[n];
		for (unsigned int j = 0; j < n; j++)
			matrix[i][j] = ref[i][j];
	}
}

FloatMatrix::~FloatMatrix()
{
	for (unsigned int i = 0; i < m; i++)
		delete matrix[i];
	delete matrix;
}


void FloatMatrix::FillFloat(float v)
{
	for (unsigned int i = 0; i < m; i++)
	{
		for (unsigned int j = 0; j < n; j++)
			matrix[i][j] = v;
	}
}

void FloatMatrix::Print()
{
	cout << "M: " << this->m << " N: " << this->n<<endl;
	for (unsigned int i = 0; i < m; i++)
	{
		for (unsigned int j = 0; j < n; j++)
			cout << fixed << matrix[i][j]<< "  ";
		cout << endl;
	}
	cout << "all" << endl;
}

void FloatMatrix::PrintRow(unsigned int a)
{
	if (a >= this->m)
		throw("Ucorrect line");

	cout << "Line: " << this->m << endl;
	for (unsigned int j = 0; j < this->n; j++)
		cout << fixed << matrix[a][j] << "  ";

	cout << endl; cout << endl;
}