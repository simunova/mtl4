#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;
    
    const unsigned n= 10;
    compressed2D<double>                           A(n, n);
    dense2D<float, mat::parameters<col_major> > B(n, n);
    morton_dense<double, 0x55555555>               C(n, n);
    morton_dense<double, 0x555555f0>               D(n, n);

    mat::hessian_setup(B, 1.0);
    mat::hessian_setup(C, 2.0);
    mat::hessian_setup(D, 3.0);
    
    std::cout << "one_norm(B) is " << one_norm(B)<< "\n";
    std::cout << "infinity_norm(B) is " << infinity_norm(B)<< "\n";
    std::cout << "frobenius_norm(B) is " << frobenius_norm(B)<< "\n";
    
    return 0;
}
