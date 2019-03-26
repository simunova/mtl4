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

#ifndef MTL_INVERT_DIAGONAL_INCLUDE
#define MTL_INVERT_DIAGONAL_INCLUDE

#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace mat {

///Returns \p A with invert diagonal
template <typename Matrix>
inline void invert_diagonal(Matrix& A)
{
    vampir_trace<3012> tracer;
    using math::reciprocal;
    namespace traits = mtl::traits;

    typename traits::row<Matrix>::type       row(A); 
    typename traits::col<Matrix>::type       col(A); 
    typename traits::value<Matrix>::type     value(A); 

    typedef typename traits::range_generator<tag::major, Matrix>::type     cursor_type;
    typedef typename traits::range_generator<tag::nz, cursor_type>::type   icursor_type;
    
    for (cursor_type cursor = begin<tag::major>(A), cend = end<tag::major>(A); cursor != cend; ++cursor) 
	for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); icursor != icend; ++icursor) 
	    if (row(*icursor) == col(*icursor))
		value(*icursor, reciprocal(value(*icursor)));
}

// Specialization for compressed2D
template <typename Value, typename Parameters>
inline void invert_diagonal(compressed2D<Value, Parameters>& A)
{
    vampir_trace<3026> tracer;
    using math::reciprocal;
    for (typename Parameters::size_type i= 0, end= std::min(num_rows(A), num_cols(A)); i < end; i++)
	A.lvalue(i, i)= reciprocal(A.lvalue(i, i));
}


}} // namespace mtl::matrix

#endif // MTL_INVERT_DIAGONAL_INCLUDE
