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

#ifndef MTL_FROBENIUS_NORM_INCLUDE
#define MTL_FROBENIUS_NORM_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/operation/max_of_sums.hpp>
#include <boost/numeric/mtl/operation/squared_abs.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace mat {

/// Frobenius norm, i.e. square root of sum of squares of all entries sqrt(sum_i sum_j(|a[i][j]|^2))
template <typename Matrix>
typename RealMagnitude<typename Collection<Matrix>::value_type>::type
inline frobenius_norm(const Matrix& matrix)
{
    vampir_trace<3010> tracer;
    using std::sqrt; using std::abs; using math::zero;
    namespace traits = mtl::traits;
    typename traits::const_value<Matrix>::type     value(matrix); 

    typedef typename Collection<Matrix>::value_type   value_type;
    typedef typename RealMagnitude<value_type>::type  real_type;
    real_type ref, sum= zero(ref);

    typedef typename traits::range_generator<tag::major, Matrix>::type     cursor_type;
    typedef typename traits::range_generator<tag::nz, cursor_type>::type   icursor_type;

    for (cursor_type cursor = begin<tag::major>(matrix), cend = end<tag::major>(matrix); cursor != cend; ++cursor) 
	for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); icursor != icend; ++icursor) 
	    sum+= squared_abs(value(*icursor));
    return sqrt(sum);
}

}} // namespace mtl::matrix

#endif // MTL_FROBENIUS_NORM_INCLUDE
