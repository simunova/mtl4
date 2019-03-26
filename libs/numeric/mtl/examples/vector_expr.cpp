// File: vector_expr.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    typedef std::complex<double>  cdouble;
    dense_vector<cdouble>         u(10), v(10);
    dense_vector<double>          w(10), x(10, 4.0);

    for (unsigned i= 0; i < size(v); i++)
	v[i]= cdouble(i+1, 10-i), w[i]= 2 * i + 2;

    u= v + w + x;
    std::cout << "u is " << u << "\n";

    u-= 3 * w;
    std::cout << "u is " << u << "\n";

    u*= 6 ;
    std::cout << "u is " << u << "\n";

    u/= 2;
    std::cout << "u is " << u << "\n";

    u+= dot(v, w) * w + 4.0 * v + 2 * w;
    std::cout << "u is " << u << "\n";

    std::cout << "i * w is " << cdouble(0,1) * w << "\n";

    return 0;
}

