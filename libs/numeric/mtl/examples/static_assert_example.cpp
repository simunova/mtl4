// Filename: static_assert_example.cpp (part of MTL4)

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

template <typename Matrix>
inline void f(const Matrix&)
{
    MTL_STATIC_ASSERT((mtl::traits::is_matrix<Matrix>::value && mtl::traits::is_dense<Matrix>::value), 
		      "f works only on dense matrices. Dork!");
    // ... now we can be sure that we have a dense matrix
}

int main()
{
    mtl::dense2D<double>      A(3, 4);
    mtl::compressed2D<double> B(3, 4);

    f(A);
    // f(B); // yields user-defined error message (try it!)

    return 0;
}
