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

#ifndef MTL_MATRIX_PERMUTATION_INCLUDE
#define MTL_MATRIX_PERMUTATION_INCLUDE

#include <boost/numeric/mtl/operation/size.hpp>
#include <boost/numeric/mtl/matrix/reorder.hpp>

namespace mtl { namespace mat {


namespace traits {

    //\ Return type of mtl::mat::permutation
    // Only for completeness	
    template <typename Value= short>
    struct permutation
    {
	typedef typename reorder<Value>::type  type;
    };
}

template <typename Value, typename PermutationVector>
typename traits::permutation<Value>::type
inline permutation(const PermutationVector& v)
{
    using mtl::size;
    return reorder(v, size(v));
}

/// Computes permutation matrix from corresponding vector
template <typename PermutationVector>
typename traits::permutation<>::type
inline permutation(const PermutationVector& v)
{
    return permutation<short>(v);
}


}} // namespace mtl::matrix

namespace mtl { namespace vec {

    /// Import into vector namespace; see \ref mtl::mat::permutation
    using mtl::mat::permutation;

}} // namespace mtl::vector

#endif // MTL_MATRIX_PERMUTATION_INCLUDE
