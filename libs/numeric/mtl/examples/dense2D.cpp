// File: dense2D.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    // A is a row-major matrix
    dense2D<double>               A(10, 10);

    // Matrices are not initialized by default
    A= 0.0;

    // Assign a value to a matrix element
    A(2, 3)= 7.0;

    // You can also use a more C-like notation
    A[2][4]= 3.0;

    std::cout << "A is \n" << A << "\n";
    
    // B is a column-major matrix
    dense2D<float, mat::parameters<tag::col_major> > B(10, 10);

    // Assign the identity matrix times 3 to B
    B= 3;
    std::cout << "B is \n" << B << "\n";

    return 0;
}

