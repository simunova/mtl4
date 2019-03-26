// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschrÃ¤nkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.
 
// #define MTL_VERBOSE_TEST

#include <iostream>
#include <boost/tuple/tuple.hpp>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/matrix/coordinate2D.hpp>

int main(int, char**)
{
    using namespace mtl;
    using mtl::io::tout;

    typedef mtl::mat::coordinate2D<double> matrix_type;
    matrix_type   B(5, 4);
    mtl::dense_vector<double> res(5, 0.0), res2(5), x(4, 1.0);
    

    tout << "num_rows = " << B.num_rows() << "\n"
	 << "num_rows = " << num_rows(B) << "\n"
	 << "num_cols = " << num_cols(B) << "\n"
	 << "size = " << size(B) << "\n"
	 << "nnz = " << nnz(B) << "\n";

    B.push_back(1, 1, 1.33);
    B.push_back(1, 2, 2.33);
    B.push_back(2, 1, 3.33);
    B.push_back(3, 1, 4.33);

    tout << "B(1, 2) = " << B(1, 2) << "\n";
    tout << "B[2][1] = " << B[2][1] << "\n";

    MTL_THROW_IF(std::abs(B(3, 1) - 4.33) > 0.001, unexpected_result());
    MTL_THROW_IF(std::abs(B(1, 2) - 2.33) > 0.001, unexpected_result());
    MTL_THROW_IF(std::abs(B(2, 1) - 2.33) > 3.001, unexpected_result());

    B.push_back(1, 0, 5.33);
    B.push_back(0, 0, 6.33);
    
    B.print_internal(tout);
    B.sort();
    
    tout << "Sorted\n";
    B.print_internal(tout);

    tout << "x=" << x << "\n" << "res=" << res << "\n";
    res = B * x ;
    tout << "res=" << res << "\n";

    res2= 6.33, 8.99, 3.33, 4.33, 0;
    res-= res;
    if (one_norm(res) > 0.1)
	throw "Matrix vector product wrong.";

    matrix_type A(5, 5, 9);
    {
	mat::inserter<matrix_type> ins(A, 3);
	ins[1][2] << 13.3;
	ins[2][2] << 23.3;
	ins[2][3] << 33.3;
	ins[2][4] << 33.3;
	ins[0][4] << 53.3;
	ins[1][4] << 6.0;
	ins[3][0] << 73.3;
    }
    tout << "A (internal) after first insertion\n";    
    A.print_internal(tout);
    MTL_THROW_IF(std::abs(A[2][3] - 33.3) > 0.001, unexpected_result());
    MTL_THROW_IF(std::abs(A[2][1]) > 0.001, unexpected_result());

    {
	mat::inserter<matrix_type, operations::update_plus<double> > ins(A);
	ins[2][3] << 3.0;
	ins[2][1] << 3.33;
    }
    tout << "A (internal) after updating insertion\n";    
    A.print_internal(tout);
    MTL_THROW_IF(std::abs(A[2][3] - 36.3) > 0.001, unexpected_result());
    MTL_THROW_IF(std::abs(A[2][1] - 3.33) > 0.001, unexpected_result());
    MTL_THROW_IF(std::abs(A[2][0]) > 0.001, unexpected_result());

    tout << "A (internal) =\n";
    A.print_internal(tout);

    tout << "A =\n" << A;
    
    traits::row<matrix_type>::type             row(A); 
    traits::col<matrix_type>::type             col(A); 
    traits::const_value<matrix_type>::type     value(A); 

    typedef traits::range_generator<tag::major, matrix_type>::type  cursor_type;

    for (cursor_type cursor = mtl::begin<tag::major>(A), cend = mtl::end<tag::major>(A); 
	 cursor != cend; ++cursor) {
	
	typedef traits::range_generator<tag::nz, cursor_type>::type icursor_type;
	for (icursor_type icursor = mtl::begin<tag::nz>(cursor), icend = mtl::end<tag::nz>(cursor); 
	     icursor != icend; ++icursor) 
	    tout << "A[" << row(*icursor) << "][" << col(*icursor) << "] = " << value(*icursor) << '\n'; 
    }

    mtl::mat::compressed2D<double> C(A);    
    tout << "C=\n"<< C << "\n";


    return 0;
}
