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

#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/glas_tag.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/maybe.hpp>
#include <boost/numeric/mtl/utility/complexity.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>
#include <boost/numeric/mtl/io/test_ostream.hpp>

using namespace std;

template <typename Matrix>
void symmetry_test(const Matrix& A)
{
    namespace traits = mtl::traits;

    typename traits::row<Matrix>::type                                 row(A); 
    typename traits::col<Matrix>::type                                 col(A); 
    typename traits::const_value<Matrix>::type                         value(A); 
    typedef typename traits::range_generator<mtl::tag::major, Matrix>::type        cursor_type;
    for (cursor_type cursor = mtl::begin<mtl::tag::major>(A), cend = mtl::end<mtl::tag::major>(A); 
	 cursor != cend; ++cursor) {
	typedef mtl::tag::nz     inner_tag;
	mtl::io::tout << "---\n";
	typedef typename traits::range_generator<inner_tag, cursor_type>::type icursor_type;
	for (icursor_type icursor = mtl::begin<inner_tag>(cursor), icend = mtl::end<inner_tag>(cursor); icursor != icend; ++icursor) {
	    mtl::io::tout << "A[" << row(*icursor) << ", " << col(*icursor) << "] = " << value(*icursor) << '\n';
	    MTL_THROW_IF(!A.indexer(A, col(*icursor), row(*icursor)), mtl::runtime_error("Symmetric counter-part not found"));    
	}
    }
} 


template <typename Matrix>
void matrix_init(Matrix& matrix)
{
    typedef typename Matrix::parameters   parameters;
    typedef typename Matrix::value_type   value_type;

    mtl::mat::compressed2D_inserter<value_type, parameters> inserter(matrix);
    inserter(1, 1) << 0;
    inserter(2, 2) << 7; inserter(1, 4) << 3; inserter(3, 2) << 9; inserter(5, 1) << 5;
}
    
 
template <typename Orientation>
void test_compressed2D(char const* name)
{
    mtl::io::tout << "\n====================\n" << name << "\n====================\n";
    typedef mtl::mat::parameters<Orientation>   parameters;
    typedef mtl::compressed2D<int, parameters>     matrix_type;
    matrix_type                                    A(6, 6); 

    matrix_init(A);
    mtl::io::tout << "A is\n" << A;
    mtl::io::tout << "Has " << A.nnz() << " non-zeros.\n";

    A.make_symmetric_pattern();
    mtl::io::tout << "A after symmetrizing is\n" << A;
    mtl::io::tout << "Has " << A.nnz() << " non-zeros.\n";
    symmetry_test(A);
}

int main(int, char**)
{
    test_compressed2D<mtl::row_major>("CRS");
    test_compressed2D<mtl::col_major>("CCS");
    
    return 0;
}
 
