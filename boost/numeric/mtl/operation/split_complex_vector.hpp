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

#ifndef MTL_VECTOR_SPLIT_COMPLEX_VECTOR_INCLUDE
#define MTL_VECTOR_SPLIT_COMPLEX_VECTOR_INCLUDE

#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/size.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace vec {

/// Split one complex-valued vector into two real-valued vectors.
/** Elements of the real vector must be assignable from the real and imaginary part of the complex elements.
    Real vectors are resized if their size is 0 otherwise the vectors must have
    the same length. **/
template <typename VectorComplex, typename VectorReal, typename VectorImaginary>
inline void split_complex_vector(const VectorComplex& c, VectorReal& r, VectorImaginary& i)
{
    using mtl::size;
    using mtl::real;
    using mtl::imag;

    vampir_trace<2034> tracer;
    r.checked_change_dim(size(c));
    i.checked_change_dim(size(c));

    for (std::size_t j= 0; j < size(c); ++j)
	r[j]= real(c[j]),
	    i[j]= imag(c[j]);
}

}} // namespace mtl::vector

#endif // MTL_VECTOR_SPLIT_COMPLEX_VECTOR_INCLUDE
