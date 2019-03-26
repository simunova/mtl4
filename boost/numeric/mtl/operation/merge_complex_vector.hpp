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

#ifndef MTL_VECTOR_MERGE_COMPLEX_VECTOR_INCLUDE
#define MTL_VECTOR_MERGE_COMPLEX_VECTOR_INCLUDE

#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace vec {

/// Merge two real-valued vectors into one complex-valued vector.
/** Elements of the complex vector must be constructible from two real elements.
    Complex vector is resized if its size is 0 otherwise the vectors must have
    the same length. **/
template <typename VectorReal, typename VectorImaginary, typename VectorComplex>
inline void merge_complex_vector(const VectorReal& r, const VectorImaginary& i, VectorComplex& c)
{
    vampir_trace<2014> tracer;
    typedef typename Collection<VectorComplex>::value_type value_type;

    MTL_THROW_IF(size(r) != size(i), incompatible_size());
    c.checked_change_dim(size(r));

    for (std::size_t j= 0; j < size(r); ++j)
	c[j]= value_type(r[j], i[j]);
}




}} // namespace mtl::vector

#endif // MTL_VECTOR_MERGE_COMPLEX_VECTOR_INCLUDE
