
#include <iostream>

#include "matrix.h"

int main()
{
    FloatMatrix tmp = FloatMatrix(2, 5);
    FloatMatrix tmp2 = FloatMatrix(2);
    tmp.FillFloat(2.0);
    tmp2.FillFloat(3.0);

    FloatMatrix tmp3 = tmp + tmp2;

    tmp3.PrintRow(5);
    cout << "Uncorrect" << endl;
    
    tmp3.Print();
   
    
}
