// File: matrix2.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

template <typename Matrix>
void fill_and_print(Matrix& A, char name)
{
    // Set values in traditional way
    A= 1.2, 3.4,
       5.6, 7.8;

    // Just print them
    std::cout << name << " is \n" << A << "\n";
}

int main(int, char**)
{
#if defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_TEMPLATE_ALIAS)
    using namespace mtl;

    // Compressed matrix
    matrix<float, compressed>                 A(2, 2);
    fill_and_print(A, 'A');
    
    // Banded matrix
    matrix<float, sparse, banded>             B(2, 2);
    fill_and_print(B, 'B');
    
    // Matrix in the ELLPACK format 
    matrix<double, ellpack>                   C(2, 2);
    fill_and_print(C, 'C');
    
    // Coordinate matrix
    matrix<float, coordinate>                 D(2, 2);
    fill_and_print(D, 'D');
    
    // Morton-order matrix with the default mask
    matrix<double, morton>                    E(2, 2);
    fill_and_print(E, 'E');
    
    // Matrix with a Morton mask is of course a Morton-order matrix
    matrix<double, mask<shark_z_64_row_mask>> F(2, 2);
    fill_and_print(F, 'F');   
#endif
    return 0;
}

