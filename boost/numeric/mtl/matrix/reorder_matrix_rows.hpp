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

#ifndef MTL_MATRIX_REORDER_MATRIX_ROWS_INCLUDE
#define MTL_MATRIX_REORDER_MATRIX_ROWS_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/operation/size.hpp>

namespace mtl { namespace mat {

/// Reorder the rows of a matrix without creating a reorder matrix
/** This is less elegant but avoids cyclic dependencies. Does not work with CCS matrices. **/
template <typename ReorderVector, typename Matrix>
Matrix reorder_matrix_rows(const ReorderVector& v, const Matrix& A)
{
    using mtl::size;

    typename mtl::traits::col<Matrix>::type                                   col(A); 
    typename mtl::traits::const_value<Matrix>::type                           value(A); 
    typedef typename mtl::traits::range_generator<tag::row, Matrix>::type     cursor_type;	
    typedef typename mtl::traits::range_generator<tag::nz, cursor_type>::type icursor_type;
    typedef typename mtl::Collection<Matrix>::size_type                       size_type;
    
    Matrix B(size(v), num_cols(A));
 
    inserter<Matrix>      ins(B, size_type(B.nnz() / num_cols(B) * 1.2));
	
    for (std::size_t i= 0; i < size(v); i++) {
	cursor_type cursor(v[i], A);   // go to row given by reorder 
	for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); icursor != icend; ++icursor) 
	    ins[i][col(*icursor)] << value(*icursor);
    }
    return B;
}

}} // namespace mtl::matrix

#endif // MTL_MATRIX_REORDER_MATRIX_ROWS_INCLUDE
