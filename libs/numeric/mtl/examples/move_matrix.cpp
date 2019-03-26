// File: move_matrix.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;

// Return a matrix with move semantics 
// (argument is only for type detection)
template <typename Matrix>
Matrix f(const Matrix&)
{
    Matrix A(3, 3);
    A= 5.0;
    return A;
}


int main(int, char**)
{
    dense2D<double>     A(3, 3), B(3, 3);
    dense2D<float>      C(3, 3);
    
    B= f(A);  // Result of f is copied shallowly, move semantic kicks in
    B= A;     // Deep copy because it's not an rvalue
    C= f(A);  // Deep copy because C has different type

    return 0;
}
