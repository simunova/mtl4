// File: vector1.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    // Define dense vector of doubles with 10 elements all set to 0.0.
    dense_vector<double>   v(10, 0.0);

    // Set element 7 to 3.0.
    v[7]= 3.0;

    std::cout << "v is " << v << "\n";
    return 0;
}
