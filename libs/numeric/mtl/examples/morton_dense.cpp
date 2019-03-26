// File: morton_dense.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    // Z-order matrix
    morton_dense<double, recursion::morton_z_mask>  A(10, 10);

    A= 0;
    A(2, 3)= 7.0;
    A[2][4]= 3.0;
    std::cout << "A is \n" << A << "\n";
    
    // B is an N-order matrix with column-major 4x4 blocks, see paper
    morton_dense<float, recursion::doppled_4_col_mask> B(10, 10);

    // Assign the identity matrix times 3 to B
    B= 3;
    std::cout << "B is \n" << B << "\n";

    return 0;
}

