// File: array_initialization.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    double array[][4]= {{3, 7.2,   0, 6}, 
			{2, 4.444, 5, 3},
			{1, 7,     9, 2}};

    dense2D<double> A(array);

    std::cout << "A = \n" << A << "\n";

    return 0;
}
