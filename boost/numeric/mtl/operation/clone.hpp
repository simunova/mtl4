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

#ifndef MTL_CLONE_INCLUDE
#define MTL_CLONE_INCLUDE

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl {

template<typename T>
struct is_clonable : boost::mpl::false_
{ };

/// Move-semantics-related anti-dot: always copy in constructor.
/** Some collections have referring semantics in copy constructors, e.g. sub-matrices.
    That means 
    \code
        Matrix B= sub_matrix(A, ...); 
    \endcode
    creates a sub-matrix of A in B. As a consequence, changes in B modify A and vice versa
    (unless it's outside the sub-matrix).
    In contrast, clone forces the copy semantics
    \code
        Matrix B= clone(sub_matrix(A, ...)); 
    \endcode
    B now contains the values of A's sub-matrix but is an independent matrix.
    Modifications to either A or B have no effect to each other.
    Requires that type T is declared clonable in terms of 
    \code
        is_clonable<T> : boost::mpl::true_ {}; 
    \endcode
**/
template <typename T>
typename boost::enable_if<is_clonable<T>, T>::type
clone(const T& x) 
{ 
    vampir_trace<3004> tracer;
    // std::cout << "Cloning clone function.\n";
    return T(x, clone_ctor()); 
}


template <typename T>
typename boost::disable_if<is_clonable<T>, T>::type
clone(const T& x) 
{ 
    // std::cout << "Not cloning clone function.\n";
    return x; 
}


} // namespace mtl

#endif // MTL_CLONE_INCLUDE
