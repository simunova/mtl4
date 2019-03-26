// Filename: ranged_for_iteration.cpp (part of MTL4)

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;
    
template <typename Matrix>
void f(Matrix& A)
{
    A= 7.0;  // Set values in diagonal
    
#if defined(MTL_WITH_AUTO) && defined(MTL_WITH_RANGEDFOR)
    // Define the property maps
    auto row=   row_map(A); 
    auto col=   col_map(A);
    auto value= const_value_map(A); 

    // Now iterate over the matrix    
    for (auto c : major_of(A))      // rows or columns
	for (auto i : nz_of(c))     // non-zeros within
	    std::cout << "A[" << row(i) << ", " << col(i) << "] = " << value(i) << '\n';    
#endif
}


int main(int, char**)
{
    // Define a row-major sparse and a column-major dense matrix
    compressed2D<double>                             A(3, 3); 
    dense2D<double, mat::parameters<col_major> >  B(3, 3); 

    f(A);
    f(B);
    
    return 0;
}
