#pragma once
#include <iostream>
#include <iomanip>
/***
Матрицы
	+Сложение/вычитание
	+Умножение на число
	+Умножение на матрицу
	+Детерминант
	-Обратная
	+Транспонированная
	+Изменение размерности
Векторы
	+Сложение/вычитание
	+Умножение на число
	+Скалярное произведение
	+-Векторное произведение
	+Длина
	+Изменение размера
***/


using namespace std;


class FloatMatrix
{
public:
	FloatMatrix();
	FloatMatrix(const FloatMatrix& obj);
	FloatMatrix(unsigned int a);
	FloatMatrix(unsigned int a, unsigned int b);//new
	FloatMatrix(unsigned int a, unsigned int b, float** ref);//new
	~FloatMatrix();

	void Print();
	void PrintRow(unsigned int a);//new
	void PrintColumn(unsigned int a);//new

	float* GetCollumn(unsigned int a);

	void FillFloat(float v);
	void SetRow(float* line, unsigned int length, unsigned int k);//new
	void SetColumn(float* column, unsigned int length, unsigned int k);//new
	void SetElement(float e, unsigned int a, unsigned int b);
	void Resize(unsigned int a, unsigned int b);
	void SwapLines(unsigned int a, unsigned int b);

	void Transpose();
	static FloatMatrix Transpose(FloatMatrix ref);
	FloatMatrix* Inverse();
	float Minor(unsigned int a, unsigned int b);//
	float Det();//new


	FloatMatrix operator + (FloatMatrix& m2)
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
				res.matrix[i][j] += m2.matrix[i][j];

		return FloatMatrix(res.m, res.n, res.matrix);

	}//new
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

	}//new
	FloatMatrix operator * (FloatMatrix& m2)
	{
		if (this->n != m2.m) {
			cout << "dfsd" << endl;
			throw "bad size";
		}
		FloatMatrix res = FloatMatrix(this->m, m2.n);
		res.FillFloat(0.0);
		for (unsigned int i = 0; i < res.m; i++)
			for (unsigned int j = 0; j < res.n; j++)
				for (unsigned int k = 0; k < this->n; k++)
					//Каждый элемент хаполняем по формуле res[i][j] = S(0,n-1)(this[i][k]*m2[k][j])
					res.matrix[i][j] += this->matrix[i][k] * m2.matrix[k][j];
		
		return FloatMatrix(res.m, res.n,res.matrix);

	};
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
	matrix = new float* [m];
	matrix[0] = new float[n];
	matrix[0][0] = 0;
}
FloatMatrix::FloatMatrix(const FloatMatrix& obj)
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

	for (unsigned int i = 0; i < this->m; i++)
		this->matrix[i][k] = column[i];
}
void FloatMatrix::SetElement(float e, unsigned int a, unsigned int b)
{
	if (a<0 || b<0 || a>this->m || b>this->n)
		throw("Uncorrect index");
	this->matrix[a][b] = e;
}

void FloatMatrix::Print()
{
	cout.precision(3);
	cout.width(9);
	cout << "M: " << this->m << " N: " << this->n << endl;
	for (unsigned int i = 0; i < this->m; i++)
	{
		for (unsigned int j = 0; j < this->n; j++)
		{
			cout << fixed << setw(9)<< this->matrix[i][j] << "  ";
		}
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
	for (unsigned int j = 0; j < this->m; j++)
		cout << fixed << matrix[j][a] << "  " << endl;

	cout << endl;
}

float* FloatMatrix::GetCollumn(unsigned int a)
{
	if (a >= this->n)
		throw("Ucorrect collumn");

	float* res = new float[this->m];

	for (unsigned int j = 0; j < this->m; j++)
		res[j] = this->matrix[j][a];

	return res;
}

void FloatMatrix::Resize(unsigned int a, unsigned int b)
{
	if (a == this->m && b == this->n)
		return;
	else
	{
		float** res = new float* [a] ;
		for (unsigned int i = 0; i < a; i++)
		{
			res[i] = new float[b];
			for (unsigned int j = 0; j < b; j++)
				if ((i < this->m) && (j < this->n))
					res[i][j] = this->matrix[i][j];
				else
					res[i][j] = 0;
		}
		this->m = a;
		this->n = b;
		delete this->matrix;
		this->matrix = res;
	}
}
void FloatMatrix::SwapLines(unsigned int a, unsigned int b)
{
	for (unsigned int i = 0; i < this->n; i++)
	{
		float tmp = this->matrix[a][i];
		this->matrix[a][i] = this->matrix[b][i];
		this->matrix[b][i] = tmp;
	}
}

void FloatMatrix::Transpose()
{
	for (unsigned int i = 0; i < this->m; i++)
		for (unsigned int j = i + 1; j < this->n; j++)
		{
			float tmp = this->matrix[i][j];
			this->matrix[i][j] = this->matrix[j][i];
			this->matrix[j][i] = tmp;
		}
			
};
FloatMatrix FloatMatrix::Transpose(FloatMatrix ref)
{
	FloatMatrix res = FloatMatrix(ref.m, ref.n);
	for (unsigned int i = 0; i < ref.m; i++)
		for (unsigned int j = i; j < ref.n; j++)
			res.matrix[i][j] = ref.matrix[j][i];

	return FloatMatrix(res.m, res.n, res.matrix);
}
//Возвращает матрицу размера (m-1;n-1) с вычеркнутый a-той строкой, b-м столбцом
float FloatMatrix::Minor(unsigned int a, unsigned int b)
{
	FloatMatrix res = FloatMatrix(this->m - 1, this->n - 1);
	bool o1 = 0, o2 = 0;
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

	return res.Det();

}
float FloatMatrix::Det()
{
	if (m != n) 
		throw("Uncorrect size");

	float det = 0.0;
	if (this->m == 1)
		return(this->matrix[0][0]);
	else if (this->m == 2)
		return((this->matrix[0][0] * this->matrix[1][1]) - (this->matrix[0][1] * this->matrix[1][0]));
	else
	{

		int k = 1;
		for (unsigned int i = 0; i < this->n; i++)
		{
			det += k * this->matrix[0][i] * this->Minor(0, i);
			k = -k;
		}
	}
	return det;


}

FloatMatrix* FloatMatrix::Inverse()
{
	if (this->Det() == 0)
		return NULL;
	FloatMatrix *res = new FloatMatrix(this->m, this->n * 2);
	res->FillFloat(0.0);
	for (unsigned int i = 0; i < this->n; i++)
	{
		float* tmp = this->GetCollumn(i);
		res->SetColumn(tmp, this->m, i);
		res->SetElement(1.0, i, i + this->n);
	}

	//Обработка первой строки
	unsigned short z = 0;
	/*while (res->matrix[0][0] == 0)
	{
		res->SwapLines(0, z);
		z++;
		if (z == this->m)
		{
			float koeff =  (res->matrix[0][0] - 1) / res->matrix[1][0];
			for (unsigned int i = 0; i < res->n; i++)
				res->matrix[0][i] -= res->matrix[1][i] * koeff;
		}
	}*/
	float koeff;
	//для всех строк
	for (unsigned int k = 0; k < this->m-1; k++) {
		for (unsigned int i = k + 1; i < res->m; i++)
		{
			if (res->matrix[k][k] == 0 && k == res->m-1)
				return NULL;
			else if (res->matrix[k][k] == 0 )
				res->SwapLines(k, k+1);
			koeff = res->matrix[i][k] / res->matrix[k][k];
			for (unsigned int j = 0; j < res->n; j++)
				res->matrix[i][j] -= res->matrix[k][j] * koeff;
		}
		
	}
	
	for (int k = this->m-1; k > 0; k--)
		for (int i = k-1; i >= 0; i--)
		{
			koeff = koeff = res->matrix[i][k] / res->matrix[k][k];
			for (unsigned int j = k; j < res->n; j++)
				res->matrix[i][j] -= res->matrix[k][j] * koeff;
		}	
	
		
	FloatMatrix* res2 = new FloatMatrix(res->m, res->m);
	for (unsigned int i = 0; i < res->m; i++)
		res2->SetColumn(res->GetCollumn(i + res->m), res->m, i);

	return res2;


}




class FloatVector
{
public:
	FloatVector();
	FloatVector(unsigned int a);
	FloatVector(unsigned int a, float e);
	FloatVector(unsigned int a, float* ref);
	FloatVector(const FloatVector& obj);


	FloatVector operator + (FloatVector& m2)
	{
		
		if (this->n != m2.n) {
			throw "Size error";
		}

		FloatVector res = FloatVector(this->n);

		for (unsigned int i = 0; i < this->n; i++)
			res.vector[i] = m2.vector[i] + this->vector[i];

		return res;

	}
	float operator*(FloatVector& m2)
	{
		if (this->n != m2.n) {
			throw "Size error";
		}
		float res = 0;
		for (unsigned int i = 0; i < this->n; i++)
			res += this->vector[i] * m2.vector[i];
		
		return res;
	}
	
	
	static FloatVector* vecMultiply(FloatVector**, unsigned short n);

	void Resize(unsigned int a);
	float Length();
	FloatVector GetOrt();

	void FillFloat(float v);
	void SetElement(float e, unsigned int a);
	void Print();

	~FloatVector();
private:
	unsigned int n;
	float* vector;

};

FloatVector::FloatVector()
{
	n = 0;
	vector = new float[0];
	for (unsigned int i = 0; i < n; i++)
		this->vector[i] = 0;
}
FloatVector::FloatVector(unsigned int a)
{
	n = a;
	vector = new float[n];
	for (unsigned int i = 0; i < n; i++)
		vector[i] = 0;
}
FloatVector::FloatVector(unsigned int a, float e) : FloatVector(a)
{
	for (unsigned int i = 0; i < this->n; i++)
		this->vector[i] = e;
}
FloatVector::FloatVector(unsigned int a, float* ref)
{
	n = a;
	vector = new float[n];
	for (unsigned int i = 0; i < n; i++)
		vector[i] = ref[i];
}
FloatVector::FloatVector(const FloatVector& obj)
{
	n = obj.n;
	vector = new float [n];
	for (unsigned int i = 0; i < n; i++)
		vector[i] = obj.vector[i];
}

void FloatVector::Resize(unsigned int a)
{
	if (a == this->n)
		return;
	else
	{
		float* newVector = new float[a];
		if (a < this->n)
			for (unsigned int i = 0; i < a; i++)
				newVector[i] = this->vector[i];
		else
		{
			for (unsigned int i = 0; i < this->n; i++)
				newVector[i] = this->vector[i];
			for (unsigned int i = this->n; i < a; i++)
				newVector[i] = 0.0;
		}

		this->n = a;
		delete this->vector;
		this->vector = newVector;
	}
		
}
float FloatVector::Length()
{
	float a = 0;
	for (unsigned int i = 0; i < this->n; i++)
		a += this->vector[i] * this->vector[i];
	return sqrt(a);
}
FloatVector FloatVector::GetOrt()
{
	float len = this->Length();
	FloatVector res = FloatVector(this->n);
	for (unsigned int i = 0; i < this->n; i++)
		res.SetElement(this->vector[i] / len, i);
	return res;
}

void FloatVector::FillFloat(float v)
{
	for (unsigned int i = 0; i < n; i++)
		vector[i] = v;
};
void FloatVector::Print()
{
	for (unsigned int i = 0; i < this->n; i++)
		cout << fixed << vector[i] << "  ";

	cout << endl; cout << endl;
};
void FloatVector::SetElement(float e, unsigned int a)
{
	vector[a] = e;
}

FloatVector* FloatVector::vecMultiply(FloatVector** vectors, unsigned short n)
{
	FloatMatrix res = FloatMatrix(n,n);
	for (unsigned int i = 0; i < n; i++)
		res.SetElement(1.0, 0, i);
	for (unsigned int i = 0; i < n-1; i++)
		res.SetRow(vectors[i]->vector, n, i + 1);
	FloatVector* ans = new FloatVector(n);
	for (unsigned int i = 0; i < n; i++)
		ans->vector[i] = res.Minor(0, i);
	return ans;
	/*FloatVector operator ^ (FloatVector& m2)
	{
		if (this->n != m2.n) {
			throw "Size error";
		}

		FloatMatrix res = FloatMatrix(3, this->n);
		for (unsigned int i = 0; i < this->n; i++)
			res.SetElement(1.0,0,i);
		res.SetRow(this->vector, this->n, 1);
		res.SetRow(m2.vector, this->n, 2);
		res.Print();

		FloatVector ans = FloatVector(this->n);
		ans.FillFloat(0.0);
		ans.Print();

		cout << res.Det();
		cout << res.Minor(0, 0) << endl;
		cout << res.Minor(0, 1) << endl;
		cout << res.Minor(0, 2) << endl;
		for (unsigned int i = 0; i < this->n; i++) {
			ans.SetElement(2.0, i);
			ans.Print();
		}
		ans.Print();
		return ans;
	}*/
}

FloatVector::~FloatVector()
{
	delete[] vector;
};



