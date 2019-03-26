// File: rich_vector_expr.cpp

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

    // Increment w by x 
    // and assign the sum of--the incremented--w and v to u
    u= v + (w+= x);
    std::cout << "u is " << u << "\n";

    // w= w * 3; x= 2; v= v + w + x; u= u + v;
    u+= v+= (w*= 3) + (x= 2);
    std::cout << "u is " << u << "w is " << w << "\n";

    return 0;
}

