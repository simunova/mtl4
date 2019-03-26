// File: recursator.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/recursion/matrix_recursator.hpp>

int main(int, char**)
{
    using namespace mtl; using std::cout;

    // Z-order matrix
    typedef morton_dense<double, recursion::morton_z_mask>  matrix_type;
    matrix_type                                             A(10, 10);
    mat::hessian_setup(A, 3.0);

    // Define a recursator over A
    mat::recursator<matrix_type>                          rec(A);

    // Access a quadrant of the matrix
    cout << "Upper right quadrant (north_east) of A is \n" << *north_east(rec) << "\n";

    // Access a quadrant's quadrant of the matrix
    cout << "Lower left (south_west) of upper right quadrant (north_east) of A is \n" 
	 << *south_west(north_east(rec)) << "\n";

    cout << "The virtual bound of 'rec' is " << rec.bound() << "\n";

    return 0;
}

