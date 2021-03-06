
#include <iostream>
#include <time.h>

#include "matrix.h"

int main()
{
    
    srand(time(0));
    unsigned int m = 4;
    float k = 1;
    float** matrix = new float*[m];
    for (unsigned int i = 0; i < m; i++)
    {
        matrix[i] = new float[m];
        for (unsigned int j = 0; j < m; j++){
            matrix[i][j] = (float)(rand()%9+1.0);
            k++;
        }
        k+=2;
            
    }
  
    Matrix<float>* tmp = new Matrix<float>(m,m, matrix);
    Matrix<float>* tmp2 = new Matrix<float>(m);
    Matrix<float>* tmp3 = new Matrix<float>(m);

    tmp->Print();
    tmp2 = tmp->Inverse();
    tmp2->Print();

    tmp3 = Matrix<float>::Multiplication(tmp, tmp2);
    tmp3->Print();

    /*FloatVector* a = new FloatVector(3);
    FloatVector* b = new FloatVector(3);
    a->FillFloat(0.0);
    b->FillFloat(0.0);
    a->SetElement(1.0, 0);
    b->SetElement(1.0, 1);

    FloatVector* c = new FloatVector(3);
    FloatVector** mass= new FloatVector* [2];
    mass[0] = a;
    mass[1] = b;

    c = FloatVector::vecMultiply(mass, 3);
    c->Print();*/


}
