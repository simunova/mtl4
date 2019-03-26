// File: vector2.cpp

#include <complex>
#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    // Define dense vector of complex with 7 elements.
    dense_vector<std::complex<float>, mtl::vec::parameters<tag::row_major> >  v(7);

    // Set all elements to 3+2i
    v= std::complex<float>(3.0, 2.0);
    std::cout << "v is " << v << "\n";

    // Set all elements to 5+0i
    v= 5.0;
    std::cout << "v is " << v << "\n";

    v= 6;
    std::cout << "v is " << v << "\n";

    return 0;
}
