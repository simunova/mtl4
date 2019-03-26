#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;
    
    const unsigned n= 10;
    dense2D<float, mat::parameters<col_major> > B(n, n);

    mat::hessian_setup(B, 1.0);
    
    std::cout << "one_norm(B) is " << one_norm(B)<< "\n";
    std::cout << "infinity_norm(B) is " << infinity_norm(B)<< "\n";
    std::cout << "frobenius_norm(B) is " << frobenius_norm(B)<< "\n";
    
    return 0;
}
