// File: recursator2.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/recursion/matrix_recursator.hpp>

int main(int, char**)
{
    using namespace mtl; using std::cout;

    typedef morton_dense<double, recursion::morton_z_mask>  matrix_type;
    matrix_type                                             A(10, 10);
    mat::hessian_setup(A, 3.0);
    mat::recursator<matrix_type>                          rec(A);

    // Create a recursator for the north_east quadrant of A
    mat::recursator<matrix_type>                          ne(north_east(rec));

    cout << "Test if recursator 'ne' refers to an empty matrix (shouldn't): " << is_empty(ne) << "\n";
    cout << "Test if north_east of 'ne' refers to an empty matrix (it should): " << is_empty(north_east(ne)) << "\n";

    cout << "Number of rows and columns of north_east quadrant is: " << num_rows(ne)
	 << " and " << num_cols(ne) << "\n";

    cout << "Test if 'ne' fills ils virtual quadrant (shouldn't): " << is_full(ne) << "\n";
    cout << "Test if north_west fills its virtual quadrant (it should): " << is_full(north_west(rec)) << "\n";

    return 0;
}

