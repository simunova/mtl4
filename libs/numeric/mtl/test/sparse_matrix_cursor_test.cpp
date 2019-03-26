// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int , char**)
{
    using namespace std;
    namespace tag= mtl::tag;
    using mtl::begin; using mtl::end; 
    
    typedef mtl::compressed2D<float>   matrix_type;
    
    matrix_type A(8, 6);
    {
	mtl::mat::inserter<matrix_type> i2(A);
	i2(2, 2) << 21; i2(2, 4) << 22; i2(6, 1) << 23; 
	i2(7, 2) << 24 << 2; i2(4, 2) << 25; i2(2, 5) << 26; 
	i2(0, 2) << 27; i2(3, 1) << 28; i2(4, 2) << 29; 
    }
    cout << "A is\n" << A << "\n";

    mtl::traits::row<matrix_type>::type             row(A); 
    mtl::traits::col<matrix_type>::type             col(A); 
    mtl::traits::const_value<matrix_type>::type     value(A); 
    typedef mtl::traits::range_generator<tag::major, matrix_type>::type  cursor_type;
    typedef mtl::traits::range_generator<tag::nz, cursor_type>::type icursor_type;

    for (cursor_type cursor = begin<tag::row>(A), cend = end<tag::row>(A); cursor != cend; ++cursor) {
	icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor);
	cout << "-----\n" << "row " << row(*icursor) << ":\n";

	for (; icursor != icend; ++icursor)
	    cout << col(*icursor) << " " << value(*icursor) << "\n";
    }
    cout << "-----\n\n";
    
    cout << "-----\ndirectly in row 2:\n";
    cursor_type cursor(2, A);
    for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); icursor != icend; ++icursor)
	cout << col(*icursor) << " " << value(*icursor) << "\n";

    cout << "A has " << A.nnz_local(1) << " non-zeros in row 1.\n";
    cout << "A has " << A.nnz_local(2) << " non-zeros in row 2.\n";

    return 0;
}
