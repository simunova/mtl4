// File: vector_norm.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    typedef std::complex<double>  cdouble;
    dense_vector<cdouble>         v(10000);

    // Initialize vector
    for (unsigned i= 0; i < size(v); i++)
	v[i]= cdouble(i+1, 10000-i);

    std::cout << "one_norm(v) is " << one_norm(v)<< "\n";
    
    std::cout << "two_norm(v) is " << two_norm(v)<< "\n";
    
    std::cout << "infinity_norm(v) is " << infinity_norm(v)<< "\n";
    
    // Unroll computation of two-norm to 6 independent statements
    std::cout << "two_norm<6>(v) is " << two_norm<6>(v)<< "\n";

    return 0;
}
