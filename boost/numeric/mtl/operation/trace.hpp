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

#ifndef MTL_TRACE_INCLUDE
#define MTL_TRACE_INCLUDE

#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace mat {

template <typename Matrix>
typename Collection<Matrix>::value_type
inline trace(const Matrix& matrix)
{
	vampir_trace<3040> tracer;
    using math::zero;
    typedef typename Collection<Matrix>::value_type value_type;

    MTL_THROW_IF(num_rows(matrix) != num_cols(matrix), matrix_not_square());

    // If matrix is empty then the result is the identity from the default-constructed value
    if (num_rows(matrix) == 0) {
	value_type ref;
	return zero(ref);
    }

    value_type value= matrix[0][0];
    for (unsigned i= 1; i < num_rows(matrix); i++)
	value+= matrix[i][i];	
    return value;
}


}} // namespace mtl::matrix

#endif // MTL_TRACE_INCLUDE
