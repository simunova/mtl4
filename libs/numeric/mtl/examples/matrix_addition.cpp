#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;
    
    const unsigned n= 10;
    compressed2D<double>                         A(n, n);
    dense2D<int, mat::parameters<col_major> > B(n, n);
    morton_dense<double, 0x555555f0>             C(n, n), D(n, n);

    mat::laplacian_setup(A, 2, 5);
    mat::hessian_setup(B, 1); mat::hessian_setup(C, 2.0); mat::hessian_setup(D, 3.0);

    D+= A - 2 * B + C;

    std::cout << "The matrices are: A=\n" << A << "B=\n" << B << "C=\n" << C << "D=\n" << D;

    return 0;
}
