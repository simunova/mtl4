// Filename: matrix_free_cg.cpp (part of MTL4)

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>


int main(int, char**)
{
    // For a more realistic example set size to 1000 or larger
    const int size = 10, N = size * size;

    typedef mtl::mat::poisson2D_dirichlet  matrix_type;
    matrix_type                               A(size, size);
    itl::pc::identity<matrix_type>            P(A);

    mtl::dense_vector<double>                 x(N, 1.0), b(N);

    b = A * x;
    x= 0;
    itl::cyclic_iteration<double>             iter(b, 100, 1.e-11, 0.0, 5);
    cg(A, x, b, P, iter);

    return 0;
}

