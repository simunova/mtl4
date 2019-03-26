// File: shallow_copy_problems_const.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;


template <typename Matrix>
void f2(Matrix& C)
{
    C= 5.0;
}

// Undermine const-ness of function argument
template <typename Matrix>
double f(const Matrix& A)
{
    Matrix B;
    B= A;
    f2(B);
    return frobenius_norm(A);
}


int main(int, char**)
{
    dense2D<double>     A(3, 3);
    A= 4.0;

    double alpha= f(A); // A is changed now!
    std::cout << alpha << "\n";

    return 0;
}
