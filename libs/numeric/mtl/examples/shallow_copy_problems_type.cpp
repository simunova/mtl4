// File: shallow_copy_problems_type.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;

int main(int, char**)
{
    dense2D<double>     A(3, 3), B(3, 3);
    dense2D<float>      C(3, 3);

    A= 4.0;

    B= A;               // Create an alias
    B*= 2.0;            // Changes also A

    C= A;               // Copies the values
    C*= 2.0;            // A is unaffected

    return 0;
}
