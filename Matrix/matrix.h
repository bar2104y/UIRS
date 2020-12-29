#pragma once
#include <iostream>
#include <iomanip>
/***
Матрицы
	+Сложение/вычитание
	+Умножение на число
	+Умножение на матрицу
	+Детерминант
	+Обратная
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

template <typename T>
class Matrix
{
public:
	// Конструкторы/деструкторы
	Matrix();
	Matrix(const Matrix& obj);
	Matrix(unsigned int a);
	Matrix(unsigned int a, unsigned int b);
	Matrix(unsigned int a, unsigned int b, T** ref);
	~Matrix();

	// Вывод информации на экран
	void Print();
	void PrintRow(unsigned int a);
	void PrintColumn(unsigned int a);

	T* GetCollumn(unsigned int a);

	// Методы добавления информации
	void FillFloat(T v);
	void SetRow(T* line, unsigned int length, unsigned int k);
	void SetColumn(T* column, unsigned int length, unsigned int k);
	void SetElement(T e, unsigned int a, unsigned int b);

	//Методы модификации матрицы
	void Resize(unsigned int a, unsigned int b);
	void SwapLines(unsigned int a, unsigned int b);

	void Transpose();
	static Matrix Transpose(Matrix ref);
	Matrix* Inverse();
	T Minor(unsigned int a, unsigned int b);
	T Det();

	//Математические операции
	void Add(Matrix<T>* m2);
	void Subtraction(Matrix<T>* m2);
	static Matrix<T>* Multiplication(Matrix<T>* m1, Matrix<T>* m2);
	void Multiplication(T k);

	/*
	Matrix operator + (Matrix& m2)
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

		Matrix res = Matrix(maxM, maxN);

		for (unsigned int i = 0; i < this->m; i++)
			for (unsigned int j = 0; j < this->n; j++)
				res.matrix[i][j] = this->matrix[i][j];

		for (unsigned int i = 0; i < m2.m; i++)
			for (unsigned int j = 0; j < m2.n; j++)
				res.matrix[i][j] += m2.matrix[i][j];

		return Matrix(res.m, res.n, res.matrix);

	}//new
	Matrix operator - (Matrix& m2)
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

		Matrix res = Matrix(maxM, maxN);

		for (unsigned int i = 0; i < this->m; i++)
			for (unsigned int j = 0; j < this->n; j++)
				res.matrix[i][j] = this->matrix[i][j];

		for (unsigned int i = 0; i < m2.m; i++)
			for (unsigned int j = 0; j < m2.n; j++)
				res.matrix[i][j] -= m2.matrix[i][j];

		return Matrix(res.m, res.n, res.matrix);

	}//new
	Matrix operator * (Matrix& m2)
	{
		if (this->n != m2.m) {
			cout << "dfsd" << endl;
			throw "bad size";
		}
		Matrix res = Matrix(this->m, m2.n);
		res.FillFloat(0.0);
		for (unsigned int i = 0; i < res.m; i++)
			for (unsigned int j = 0; j < res.n; j++)
				for (unsigned int k = 0; k < this->n; k++)
					//Каждый элемент хаполняем по формуле res[i][j] = S(0,n-1)(this[i][k]*m2[k][j])
					res.matrix[i][j] += this->matrix[i][k] * m2.matrix[k][j];
		
		return Matrix(res.m, res.n,res.matrix);

	};*/
	Matrix operator * (int k)
	{
		for (unsigned int i = 0; i < this->m; i++)
			for (unsigned int j = 0; j < this->n; j++)
				this->matrix[i][j] *= k;
		return Matrix(this->m, this->n, this->matrix);
	}
	Matrix operator * (T k)
	{
		for (unsigned int i = 0; i < this->m; i++)
			for (unsigned int j = 0; j < this->n; j++)
				this->matrix[i][j] *= k;
		return Matrix(this->m, this->n, this->matrix);
	}
	Matrix operator * (double k)
	{
		for (unsigned int i = 0; i < this->m; i++)
			for (unsigned int j = 0; j < this->n; j++)
				this->matrix[i][j] *= T(k);
		return Matrix(this->m, this->n, this->matrix);
	}

private:
	unsigned int m, n;//строка, столбец
	T** matrix;
};

/// <summary>
/// Конструктор по умолчанию
/// </summary>
/// <typeparam name="T">Базовый тип ячеек матрицы</typeparam>
template <typename T>
Matrix<T>::Matrix()
{
	n = 1;
	m = 1;
	matrix = new T* [m];
	matrix[0] = new T[n];
	matrix[0][0] = 0;
}
/// <summary>
/// Глубокий конструктор
/// </summary>
/// <param name="obj">Экземпляр существующего объекта</param>
template <typename T>
Matrix<T>::Matrix(const Matrix& obj)
{
	m = obj.m;
	n = obj.n;
	matrix = new T* [m];
	for (unsigned int i = 0; i < m; i++)
	{
		matrix[i] = new T[n];
		for (unsigned int j = 0; j < n; j++)
			matrix[i][j] = 0;
	}

}
/// <summary>
/// Конструктор квадратной матрицы, заполненной нулями
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="a">Размерность матрицы</param>
template <typename T>
Matrix<T>::Matrix(unsigned int a)
{
	n = a;
	m = a;
	matrix = new T* [m];
	for (unsigned int i = 0; i < m; i++)
	{
		matrix[i] = new T[n];
		for (unsigned int j = 0; j < n; j++)
			matrix[i][j] = 0;
	}
}
/// <summary>
/// Конструктор матрицы, заполненной нулями
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="a">Количество строк</param>
/// <param name="b">Количество столбцов</param>
template <typename T>
Matrix<T>::Matrix(unsigned int a, unsigned int b)
{
	m = a;
	n = b;
	matrix = new T* [m];
	for (unsigned int i = 0; i < m; i++)
	{
		matrix[i] = new T[n];
		for (unsigned int j = 0; j < n; j++)
			matrix[i][j] = 0;
	}
}
/// <summary>
/// Создание матрицы с референсного массива
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="a">Количество строк</param>
/// <param name="b">Количество столбцов</param>
/// <param name="ref">Референсная матрица</param>
template <typename T>
Matrix<T>::Matrix(unsigned int a, unsigned int b, T** ref)
{
	m = a;
	n = b;
	matrix = new T* [m];
	for (unsigned int i = 0; i < m; i++)
	{
		matrix[i] = new T[n];
		for (unsigned int j = 0; j < n; j++)
			matrix[i][j] = ref[i][j];
	}
}
template <typename T>
Matrix<T>::~Matrix()
{
	for (unsigned int i = 0; i < m; i++)
		delete matrix[i];
	delete matrix;
}



template <typename T>
void Matrix<T>::FillFloat(T v)
{
	for (unsigned int i = 0; i < m; i++)
	{
		for (unsigned int j = 0; j < n; j++)
			matrix[i][j] = v;
	}
}
template <typename T>
void Matrix<T>::SetRow(T* line, unsigned int length, unsigned int k)
{
	if (length != this->n)
		throw("Uncorrect length");

	for (unsigned int i = 0; i < n; i++)
		this->matrix[k][i] = line[i];
}
template <typename T>
void Matrix<T>::SetColumn(T* column, unsigned int length, unsigned int k)
{
	if (length != this->m)
		throw("Uncorrect length");

	for (unsigned int i = 0; i < this->m; i++)
		this->matrix[i][k] = column[i];
}
template <typename T>
void Matrix<T>::SetElement(T e, unsigned int a, unsigned int b)
{
	if (a<0 || b<0 || a>this->m || b>this->n)
		throw("Uncorrect index");
	this->matrix[a][b] = e;
}


template <typename T>
void Matrix<T>::Print()
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
template <typename T>
void Matrix<T>::PrintRow(unsigned int a)
{
	if (a >= this->m)
		throw("Ucorrect line");

	cout << "Line: " << a << endl;
	for (unsigned int j = 0; j < this->n; j++)
		cout << fixed << matrix[a][j] << "  ";

	cout << endl; cout << endl;
}
template <typename T>
void Matrix<T>::PrintColumn(unsigned int a)
{
	if (a >= this->n)
		throw("Ucorrect line");

	cout << "Column: " << a << endl;
	for (unsigned int j = 0; j < this->m; j++)
		cout << fixed << matrix[j][a] << "  " << endl;

	cout << endl;
}

template <typename T>
T* Matrix<T>::GetCollumn(unsigned int a)
{
	if (a >= this->n)
		throw("Ucorrect collumn");

	T* res = new T[this->m];

	for (unsigned int j = 0; j < this->m; j++)
		res[j] = this->matrix[j][a];

	return res;
}
template <typename T>
void Matrix<T>::Resize(unsigned int a, unsigned int b)
{
	if (a == this->m && b == this->n)
		return;
	else
	{
		T** res = new T* [a] ;
		for (unsigned int i = 0; i < a; i++)
		{
			res[i] = new T[b];
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
template <typename T>
void Matrix<T>::SwapLines(unsigned int a, unsigned int b)
{
	for (unsigned int i = 0; i < this->n; i++)
	{
		T tmp = this->matrix[a][i];
		this->matrix[a][i] = this->matrix[b][i];
		this->matrix[b][i] = tmp;
	}
}

template <typename T>
void Matrix<T>::Transpose()
{
	for (unsigned int i = 0; i < this->m; i++)
		for (unsigned int j = i + 1; j < this->n; j++)
		{
			T tmp = this->matrix[i][j];
			this->matrix[i][j] = this->matrix[j][i];
			this->matrix[j][i] = tmp;
		}
			
};
template <typename T>
Matrix<T> Matrix<T>::Transpose(Matrix ref)
{
	Matrix res = Matrix(ref.m, ref.n);
	for (unsigned int i = 0; i < ref.m; i++)
		for (unsigned int j = i; j < ref.n; j++)
			res.matrix[i][j] = ref.matrix[j][i];

	return Matrix(res.m, res.n, res.matrix);
}
//Возвращает матрицу размера (m-1;n-1) с вычеркнутый a-той строкой, b-м столбцом
template <typename T>
T Matrix<T>::Minor(unsigned int a, unsigned int b)
{
	Matrix res = Matrix(this->m - 1, this->n - 1);
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
template <typename T>
T Matrix<T>::Det()
{
	if (m != n) 
		throw("Uncorrect size");

	T det = 0.0;
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
template <typename T>
Matrix<T>* Matrix<T>::Inverse()
{
	if (this->Det() == 0)
		return NULL;
	Matrix *res = new Matrix(this->m, this->n * 2);
	res->FillFloat(0.0);
	for (unsigned int i = 0; i < this->n; i++)
	{
		T* tmp = this->GetCollumn(i);
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
			T koeff =  (res->matrix[0][0] - 1) / res->matrix[1][0];
			for (unsigned int i = 0; i < res->n; i++)
				res->matrix[0][i] -= res->matrix[1][i] * koeff;
		}
	}*/
	T koeff;
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
	
	for (unsigned int k = this->m-1; k > 0; k--)
		for (int i = k-1; i >= 0; i--)
		{
			koeff = koeff = res->matrix[i][k] / res->matrix[k][k];
			for (unsigned int j = k; j < res->n; j++)
				res->matrix[i][j] -= res->matrix[k][j] * koeff;
		}

	for (unsigned int i = 0; i < this->m;  i++)
	{
		koeff =  res->matrix[i][i];
		for (unsigned int j = 0; j < res->n; j++)
			res->matrix[i][j] = res->matrix[i][j] / koeff;
	}
	
		
	Matrix* res2 = new Matrix(res->m, res->m);
	for (unsigned int i = 0; i < res->m; i++)
		res2->SetColumn(res->GetCollumn(i + res->m), res->m, i);

	return res2;


}



template <typename T>
void Matrix<T>::Add(Matrix<T>* m2)
{
	unsigned int maxM, maxN, minM, minN;
	if (this->m != m2.m or this->n != m2.m)
	{
		cout << "Неправильная размерность при сложении" << endl;
		return;
	}

	for (unsigned int i = 0; i < this.m; i++)
		for (unsigned int j = 0; j < this.n; j++)
			this.matrix[i][j] += m2.matrix[i][j];

	return;
}
template <typename T>
void Matrix<T>::Subtraction(Matrix<T>* m2)
{
	unsigned int maxM, maxN, minM, minN;
	if (this->m != m2.m or this->n != m2.m)
	{
		cout << "Неправильная размерность при сложении" << endl;
		return;
	}

	for (unsigned int i = 0; i < this.m; i++)
		for (unsigned int j = 0; j < this.n; j++)
			this.matrix[i][j] -= m2.matrix[i][j];

	return;
}
template <typename T>
Matrix<T>* Matrix<T>::Multiplication(Matrix<T>* m1, Matrix<T>* m2)
{
	if (m1->n != m2->m) {
		cout << "Неправильная размерность при умножении" << endl;
		return NULL;
	}
	Matrix<T>* res = new Matrix<T>(m1->m, m2->n);
	res->FillFloat(0.0);
	for (unsigned int i = 0; i < res->m; i++)
		for (unsigned int j = 0; j < res->n; j++)
			for (unsigned int k = 0; k < m1->n; k++)
				//Каждый элемент хаполняем по формуле res[i][j] = S(0,n-1)(this[i][k]*m2[k][j])
				res->matrix[i][j] += m1->matrix[i][k] * m2->matrix[k][j];
	

	return res;
}
template <typename T>
void Matrix<T>::Multiplication(T k)
{
	for (unsigned int i = 0; i < this->m; i++)
		for (unsigned int j = 0; j < this->n; j++)
			this->matrix[i][j] *= k;
	return;
}






/**********************************************************************************************/
template <typename T>
class Vector
{
public:
	Vector();
	Vector(unsigned int a);
	Vector(unsigned int a, T e);
	Vector(unsigned int a, T* ref);
	Vector(const Vector<T>& obj);


	Vector operator + (Vector& m2)
	{
		
		if (this->n != m2.n) {
			throw "Size error";
		}

		Vector res = Vector(this->n);

		for (unsigned int i = 0; i < this->n; i++)
			res.vector[i] = m2.vector[i] + this->vector[i];
		
		return res;

	}
	T operator*(Vector& m2)
	{
		if (this->n != m2.n) {
			throw "Size error";
		}
		T res = 0;
		for (unsigned int i = 0; i < this->n; i++)
			res += this->vector[i] * m2.vector[i];
		
		return res;
	}
	
	Matrix<T> ToMatrix()
	{
		Matrix<T>* res = new Matrix<T>(1, this->n);
		for (unsigned int i = 0; i < this->n; i++)
			res->matrix[1][i] = this->vector[i];
	};
	static Vector* vecMultiply(Vector**, unsigned short n);

	void Resize(unsigned int a);
	T Length();
	Vector GetOrt();

	void FillFloat(T v);
	void SetElement(T e, unsigned int a);
	void Print();

	~Vector();
private:
	unsigned int n;
	T* vector;

};
template <typename T>
Vector<T>::Vector()
{
	n = 0;
	vector = new T[0];
	for (unsigned int i = 0; i < n; i++)
		this->vector[i] = 0;
}
template <typename T>
Vector<T>::Vector(unsigned int a)
{
	n = a;
	vector = new T[n];
	for (unsigned int i = 0; i < n; i++)
		vector[i] = 0;
}
template <typename T>
Vector<T>::Vector(unsigned int a, T e) : Vector(a)
{
	for (unsigned int i = 0; i < this->n; i++)
		this->vector[i] = e;
}
template <typename T>
Vector<T>::Vector(unsigned int a, T* ref)
{
	n = a;
	vector = new T[n];
	for (unsigned int i = 0; i < n; i++)
		vector[i] = ref[i];
}
template <typename T>
Vector<T>::Vector(const Vector& obj)
{
	n = obj.n;
	vector = new T [n];
	for (unsigned int i = 0; i < n; i++)
		vector[i] = obj.vector[i];
}

template <typename T>
void Vector<T>::Resize(unsigned int a)
{
	if (a == this->n)
		return;
	else
	{
		T* newVector = new T[a];
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
template <typename T>
T Vector<T>::Length()
{
	T a = 0;
	for (unsigned int i = 0; i < this->n; i++)
		a += this->vector[i] * this->vector[i];
	return sqrt(a);
}
template <typename T>
Vector<T> Vector<T>::GetOrt()
{
	T len = this->Length();
	Vector res = Vector(this->n);
	for (unsigned int i = 0; i < this->n; i++)
		res.SetElement(this->vector[i] / len, i);
	return res;
}

template <typename T>
void Vector<T>::FillFloat(T v)
{
	for (unsigned int i = 0; i < n; i++)
		vector[i] = v;
}
template <typename T>
void Vector<T>::Print()
{
	for (unsigned int i = 0; i < this->n; i++)
		cout << fixed << vector[i] << "  ";

	cout << endl; cout << endl;
}
template <typename T>
void Vector<T>::SetElement(T e, unsigned int a)
{
	vector[a] = e;
}

template <typename T>
Vector<T>* Vector<T>::vecMultiply(Vector<T>** vectors, unsigned short n)
{
	Matrix<T> res = Matrix<T>(n,n);
	for (unsigned int i = 0; i < n; i++)
		res.SetElement(1.0, 0, i);
	for (unsigned int i = 0; i < n-1; i++)
		res.SetRow(vectors[i]->vector, n, i + 1);
	Vector* ans = new Vector(n);
	for (unsigned int i = 0; i < n; i++)
		ans->vector[i] = res.Minor(0, i);
	return ans;
	/*Vector operator ^ (Vector& m2)
	{
		if (this->n != m2.n) {
			throw "Size error";
		}

		Matrix res = Matrix(3, this->n);
		for (unsigned int i = 0; i < this->n; i++)
			res.SetElement(1.0,0,i);
		res.SetRow(this->vector, this->n, 1);
		res.SetRow(m2.vector, this->n, 2);
		res.Print();

		Vector ans = Vector(this->n);
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

template <typename T>
Vector<T>::~Vector()
{
	delete[] vector;
};



