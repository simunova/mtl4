#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;
    
    const unsigned n= 10;
    dense2D<int, mat::parameters<col_major> > b(n, n);
    morton_dense<double, 0x55555555>             c(n, n);
    morton_dense<double, 0x555555f0>             d(n, n);

    b= 0; d= 0;

    d(2, 2)= 3; 
    b[2][3]= 4;

    c= 2.0;

    std::cout << std::complex<double>(0, 1) * c;
    d*= 7.0;

    d+= c - b;
    d-= b;

    return 0;
}
