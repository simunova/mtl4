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

#ifndef MTL_MATRIX_REORDER_REF_INCLUDE
#define MTL_MATRIX_REORDER_REF_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/size.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>

namespace mtl { namespace mat {


template <typename ReorderVector, typename Matrix>
void reorder_ref(const ReorderVector& v, Matrix& A, std::size_t cols= 0)
{
    using math::one; using mtl::size;
    typedef typename Collection<Matrix>::value_type value_type;

    if (size(v) == 0) {
	A.change_dim(0, cols); return; }

    // Find maximal entry (don't use mtl::max to allow for arrays and others)
    std::size_t  s= size(v),
	         my_max= std::size_t(*std::max_element(&v[0], &v[0] + s)) + 1;

    if (cols == 0) 
	cols= my_max;
    else
	MTL_THROW_IF(my_max > cols, range_error("Too large value in reorder vector"));

    A.change_dim(s, cols);
    inserter<Matrix>      ins(A, 1);
    for (std::size_t i= 0; i < s; i++)
	ins[i][v[i]] << one(value_type());
}


}} // namespace mtl::matrix

#endif // MTL_MATRIX_REORDER_REF_INCLUDE
