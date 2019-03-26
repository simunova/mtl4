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

#ifndef MTL_INDEX_INCLUDE
#define MTL_INDEX_INCLUDE

// The whole idea of changing indices is insane!
// Thus, the file shouldn't exist at all.

#include <boost/mpl/if.hpp>

namespace mtl { namespace index {

// Index like in C (identical with internal representation)
struct c_index {};

// Index like Fortran
struct f_index {};


#if 0
// Which index has type T
template <class T> struct which_index
{
    typedef typename boost::mpl::if_c<
          traits::is_mtl_type<T>::value
        , typename T::index_type   // mtl data shall know their type
        , c_index                  // others are by default c
        >::type type;
};
#endif


template <class T> struct which_index
{
    typedef typename T::index_type type;
};

// Change from internal representation to requested index type
template <class T> inline T change_to(c_index, T i) 
{
    return i; 
}

template <class T> inline T change_to(f_index, T i) 
{ 
    return i + 1; 
}

// Change from requested index type to internal representation
template <class T> inline T change_from(c_index, T i) 
{ 
    return i; 
}

template <class T> inline T change_from(f_index, T i) 
{ 
    return i - 1; 
}

}} // namespace mtl::index

#endif // MTL_INDEX_INCLUDE
