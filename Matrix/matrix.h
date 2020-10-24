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


	void Print();
	void PrintRow(unsigned int a);
	void PrintColumn(unsigned int a);


	void FillFloat(float v);
	void SetRow(float* line, unsigned int length, unsigned int k);
	void SetColumn(float* column, unsigned int length, unsigned int k);
	void SetElement(float e, unsigned int a, unsigned int b);


	float Minor(unsigned int a, unsigned int b);
	float Det();


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
			for (unsigned int j = 0; j < this->n; j++)
				res.matrix[i][j] = this->matrix[i][j];

		for (unsigned int i = 0; i < m2.m; i++)
			for (unsigned int j = 0; j < m2.n; j++)
				res.matrix[i][j] += m2.matrix[i][j];
		
		return FloatMatrix(res.m,res.n, res.matrix);

	}
	FloatMatrix operator - (FloatMatrix& m2)
	{
		unsigned int maxM, maxN, minM, minN;
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
			for (unsigned int j = 0; j < this->n; j++)
				res.matrix[i][j] = this->matrix[i][j];

		for (unsigned int i = 0; i < m2.m; i++)
			for (unsigned int j = 0; j < m2.n; j++)
				res.matrix[i][j] -= m2.matrix[i][j];

		return FloatMatrix(res.m, res.n, res.matrix);

	}
	FloatMatrix operator * (int k)
	{
		for (unsigned int i = 0; i < this->m; i++)
			for (unsigned int j = 0; j < this->n; j++)
				this->matrix[i][j] *= k;
		return FloatMatrix(this->m, this->n, this->matrix);
	}
	FloatMatrix operator * (float k)
	{
		for (unsigned int i = 0; i < this->m; i++)
			for (unsigned int j = 0; j < this->n; j++)
				this->matrix[i][j] *= k;
		return FloatMatrix(this->m, this->n, this->matrix);
	}
	FloatMatrix operator * (double k)
	{
		for (unsigned int i = 0; i < this->m; i++)
			for (unsigned int j = 0; j < this->n; j++)
				this->matrix[i][j] *= float(k);
		return FloatMatrix(this->m, this->n, this->matrix);
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

void FloatMatrix::SetRow(float* line, unsigned int length, unsigned int k)
{
	if (length != this->n)
		throw("Uncorrect length");

	for (unsigned int i = 0; i < n; i++)
		this->matrix[k][i] = line[i];
}

void FloatMatrix::SetColumn(float* column, unsigned int length, unsigned int k)
{
	if (length != this->m)
		throw("Uncorrect length");

	for (unsigned int i = 0; i < n; i++)
		this->matrix[i][k] = column[i];
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

	cout << "Line: " << a << endl;
	for (unsigned int j = 0; j < this->n; j++)
		cout << fixed << matrix[a][j] << "  ";

	cout << endl; cout << endl;
}

void FloatMatrix::PrintColumn(unsigned int a)
{
	if (a >= this->n)
		throw("Ucorrect line");

	cout << "Column: " << a << endl;
	for (unsigned int j = 0; j < this->n; j++)
		cout << fixed << matrix[j][a] << "  " << endl;

	cout << endl; 
}

void FloatMatrix::SetElement(float e, unsigned int a, unsigned int b)
{
	if (a<0 || b<0 || a>this->m || b>this->n)
		throw("Uncorrect index");
	this->matrix[a][b] = e;
}


//Возвращает матрицу размера (m-1;n-1) с вычеркнутый a-той строкой, b-м столбцом
float FloatMatrix::Minor(unsigned int a, unsigned int b)
{
	FloatMatrix res = FloatMatrix(this->m-1, this->n-1);
	bool o1=0, o2=0;
	for (unsigned int i = 0; i < m; i++)
	{
		if (i == a)
		{
			o1 = 1;
			continue;
		}
		o2 = 0;
		for (unsigned int j = 0; j < n; j++)
		{
			if (j == b)
			{
				o2 = true;
				continue;
			}
			res.matrix[i - o1][j - o2] = this->matrix[i][j];
		}
	}
	cout << "DETOPR   " << res.Det() << endl;
	return res.Det();

}

float FloatMatrix::Det()
{
	if (m != n)
		throw("Uncorrect size");
	float det = 0;
	if (this->m == 1)
		return(this->matrix[0][0]);
	else if (this->m == 2)
		return((this->matrix[0][0] * this->matrix[1][1]) - (this->matrix[0][1] * this->matrix[1][0]));
	else
	{
		
		int k = 1;
		for (unsigned int i = 0; i < this->n; i++)
		{
			det += k*this->matrix[0][i] * this->Minor(0, i);
			k = -k;
		}
	}
	return det;
	
	
}