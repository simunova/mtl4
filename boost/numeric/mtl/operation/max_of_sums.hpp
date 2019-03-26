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

#ifndef MTL_MAX_OF_SUMS_INCLUDE
#define MTL_MAX_OF_SUMS_INCLUDE

#include <boost/numeric/mtl/concept/magnitude.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/linear_algebra/operators.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#include <numeric>
#include <cmath>


namespace mtl { namespace impl {

// We need property map of the minor index 
template <typename Matrix, typename MinorIndex>
typename RealMagnitude<typename Collection<Matrix>::value_type>::type
inline max_of_sums(const Matrix& matrix, bool aligned, MinorIndex minor_index, std::size_t dim2)
{
    vampir_trace<2012> tracer;
    using std::max; using std::abs; using math::zero;

    typedef typename Collection<Matrix>::value_type   value_type;
    typedef typename RealMagnitude<value_type>::type  real_type;
    real_type ref, my_zero= zero(ref);

    // If matrix is empty then the result is the identity from the default-constructed value
    if (num_rows(matrix) == 0 || num_cols(matrix) == 0)
	return my_zero;

    typedef typename traits::range_generator<tag::major, Matrix>::type     cursor_type;
    typedef typename traits::range_generator<tag::nz, cursor_type>::type   icursor_type;
    typename traits::const_value<Matrix>::type                             value(matrix); 

    if (aligned) {
	real_type maxv= my_zero;
	for (cursor_type cursor = begin<tag::major>(matrix), cend = end<tag::major>(matrix); cursor != cend; ++cursor) {
	    real_type sum= my_zero;
	    for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); icursor != icend; ++icursor)
		sum+= abs(value(*icursor));
	    maxv= max(maxv, sum);
	}
	return maxv;
    }

    // If matrix has other orientation, we compute all sums in a vector
    dense_vector<real_type>   sums(dim2, my_zero);
    for (cursor_type cursor = begin<tag::major>(matrix), cend = end<tag::major>(matrix); cursor != cend; ++cursor)
	for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); icursor != icend; ++icursor)
	    sums[minor_index(*icursor)]+= abs(value(*icursor));
    // replace by mtl::accumulate<8>
    return std::accumulate(sums.begin(), sums.end(), my_zero, math::max<real_type>());
}


}} // namespace mtl::impl

#endif // MTL_MAX_OF_SUMS_INCLUDE
