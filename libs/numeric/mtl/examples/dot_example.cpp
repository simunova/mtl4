// Filename: dot_example.cpp (part of MTL4)

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    typedef std::complex<double>  cdouble;
    dense_vector<cdouble>         v(10000), x(10, cdouble(3, 2));
    dense_vector<double>          w(10000);

    for (unsigned i= 0; i < size(v); i++)
	v[i]= cdouble(i+1, 10000-i), w[i]= 2 * i + 2;

    std::cout << "dot(v, w) is " << dot(v, w)<< "\n";
    
    std::cout << "dot<6>(v, w) is " <<  dot<6>(v, w)<< "\n";

    std::cout << "dot_real<6>(v, w) is " <<  dot_real<6>(v, w)<< "\n";
    
    std::cout << "conj(x) is " <<  mtl::conj(x)<< "\n"; // ADL doesn't work here in g++ 4.4

    return 0;
}

