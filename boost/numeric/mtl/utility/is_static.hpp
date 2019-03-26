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

#ifndef MTL_TRAITS_IS_STATIC_INCLUDE
#define MTL_TRAITS_IS_STATIC_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace traits {

    /// Meta-function whether a certain type has static size
    template <typename T> struct is_static : boost::mpl::false_ {};

    template <std::size_t Size> struct is_static<mtl::vec::fixed::dimension<Size> > : boost::mpl::true_ {};
    template <std::size_t Rows, std::size_t Cols> struct is_static<mtl::fixed::dimensions<Rows, Cols> > : boost::mpl::true_ {};

    template <typename V, typename P> struct is_static<mtl::vec::dense_vector<V, P> > : is_static<typename P::dimension> {};

    template <typename V, typename P> struct is_static<mtl::mat::dense2D<V, P> > : is_static<typename P::dimensions> {};
    template <typename V, std::size_t M, typename P> struct is_static<mtl::mat::morton_dense<V, M, P> > : is_static<typename P::dimensions> {};
    template <typename V, typename P> struct is_static<mtl::mat::compressed2D<V, P> > : is_static<typename P::dimensions> {};

}} // namespace mtl::traits

#endif // MTL_TRAITS_IS_STATIC_INCLUDE
