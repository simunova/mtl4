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

#ifndef MTL_VECTOR_LAZY_REDUCTION_INCLUDE
#define MTL_VECTOR_LAZY_REDUCTION_INCLUDE

namespace mtl { namespace vec {

/// Helper class for lazy evaluation
template <typename Vector, typename Functor>
struct lazy_reduction
{
    lazy_reduction(const Vector& v) : v(v) {}
    const Vector& v;
};

template <typename Vector, typename Functor>
inline std::size_t size(const lazy_reduction<Vector, Functor>& lazy)
{   
    return size(lazy.v); 
}

}} // namespace mtl::vector

#endif // MTL_VECTOR_LAZY_REDUCTION_INCLUDE
