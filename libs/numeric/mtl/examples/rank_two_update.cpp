#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;
    typedef std::complex<double>  cdouble;
    
    const unsigned n= 8;
    dense2D<cdouble>              A(n, n);

    A= 3.0;

    dense_vector<cdouble>         v(n), w(n);
    for (unsigned i= 0; i < size(v); i++)
	v[i]= cdouble(i+1, n-i), w[i]= cdouble(i+n);

    rank_one_update(A, v, w);
    std::cout << "A after rank-one update is \n" 
	      << with_format(A, 9, 3) << "\n";

    A= 3.0;
    rank_two_update(A, v, w);
    std::cout << "A after rank-two update is \n"
	      << with_format(A, 9, 3) << "\n";

    return 0;
}
