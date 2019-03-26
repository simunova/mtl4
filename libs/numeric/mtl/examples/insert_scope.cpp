// Filename: insert_scope.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;

int main(int, char**)
{
    compressed2D<double>              A(3, 3);
    {
	mat::inserter<compressed2D<double> > ins(A);  
	ins[0][0] << 2.0;
	ins[1][2] << 0.5;
	ins[2][1] << 3.0;
    } // ins is destroyed here

    std::cout << "A is \n" << A;  // we can access A now

    return 0;
}

