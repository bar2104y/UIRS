
#include <iostream>

#include "matrix.h"

int main()
{
    
    //FloatMatrix tmp = FloatMatrix(2, 5);
    int k = 1;
    int n = 3;
    float row [3] = {1.0,2.0,3.0};
    float** matrix = new float* [n];
    for (unsigned int i = 0; i < n; i++)
    {
        matrix[i] = new float[n];
        for (unsigned int j = 0; j < n; j++) {
            matrix[i][j] = k;
            k++;
        }
    }
    FloatMatrix tmp = FloatMatrix(n, n, matrix);
    tmp.FillFloat(2.0);
    tmp.Print();
    tmp.SetRow(row, 3, 1);
    tmp.Print();
    tmp.SetColumn(row, 3, 1);
    tmp.SetElement(8.0, 1, 1);
    tmp.Print();


    float d = tmp.Det();
    cout << d << endl;

    /*FloatMatrix tmp = FloatMatrix(3);
    FloatMatrix tmp2 = FloatMatrix(5);
    
    tmp.FillFloat(2.0);
    tmp2.FillFloat(3.0);

    FloatMatrix tmp3 = tmp2 * 2;
    FloatMatrix tmp4 = (tmp2 + tmp)*2;
    tmp.Print();
    tmp2.Print();
    tmp3.Print();
    tmp4.Print();*/

   
    
}
