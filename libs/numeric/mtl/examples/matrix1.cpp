// File: matrix1.cpp

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

    // Create a 2 by 2 dense matrix
    matrix<double>                            A(2, 2);
    fill_and_print(A, 'A');
    
    // Now a sparse matrix
    matrix<double, sparse>                    B(2, 2);
    fill_and_print(B, 'B');
    
    // A column-major sparse matrix (CCS by default)
    matrix<double, sparse, column_major>      C(2, 2);
    fill_and_print(C, 'C');
    
    // A Fortran-like matrix
    matrix<double, col_major>                 D(2, 2);
    fill_and_print(D, 'D');
    
    // A sparse matrix with a shorter size type
    matrix<float, sparse, as_size_type<int> > E(2, 2);
    fill_and_print(E, 'E');
    
    // A (dense) matrix with compile-time size
    matrix<double, dim<2, 2> >                F;
    fill_and_print(F, 'F');
#endif
    return 0;
}

