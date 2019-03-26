// File: vector_min_max.cpp

#include <iostream>
#include <cmath>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using mtl::max; using std::pow; // to avoid possible ambiguity with ARPREC

    mtl::dense_vector<double>         v(100);

    for (unsigned i= 0; i < size(v); i++)
	v[i]= double(i+1) * pow(-1.0, int(i)); // Amb. in MSVC

    std::cout << "max(v) is " << max(v)<< "\n";
    
    std::cout << "min(v) is " <<  min(v)<< "\n";
    
    std::cout << "max<6>(v) is " <<  max<6>(v)<< "\n";

    return 0;
}

