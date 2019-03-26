// Filename: complete_iteration.cpp (part of MTL4)

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;
    
template <typename Matrix>
void f(Matrix& A)
{
    using traits::range_generator; 
    A= 7.0;    // Set values in diagonal
 
    // Types of outer and inner cursor
    typedef typename range_generator<tag::major, Matrix>::type c_type;
    typedef typename range_generator<tag::nz, c_type>::type    ic_type;

    // Define the property maps
    typename traits::row<Matrix>::type               row(A); 
    typename traits::col<Matrix>::type               col(A);
    typename traits::const_value<Matrix>::type       value(A); 

    // Now iterate over the matrix    
    for (c_type cursor= begin<tag::major>(A), cend= end<tag::major>(A); cursor != cend; ++cursor)
       for (ic_type icursor= begin<tag::nz>(cursor), icend= end<tag::nz>(cursor); icursor != icend; ++icursor)
	   std::cout << "A[" << row(*icursor) << ", " << col(*icursor) << "] = " << value(*icursor) << '\n';    
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
