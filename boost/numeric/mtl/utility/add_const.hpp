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

#ifndef MTL_TRAITS_ADD_CONST_INCLUDE
#define MTL_TRAITS_ADD_CONST_INCLUDE

namespace mtl { namespace traits {

/// Add const to data
/** In case of pointers, constify data of first level (even if it is a pointer itself)
    instead of (outer) address.
    Other types are constified like in boost::add_const.
    \sa add_const_to_root, add_const_to_all  **/
template <typename T>
struct add_const_to_data
{
    typedef T const  type;
};

template <typename T*>
struct add_const_to_data
{
    typedef T const *  type;
};

/// Add const to data at the root
/** In case of pointers, constify data of innermost level and leave other levels as they are.
    instead of (outer) address.
    Other types are constified like in boost::add_const.
    \sa add_const_to_data, add_const_to_all  **/
template <typename T>
struct add_const_to_root
{
    typedef T const  type;
};

template <typename T*>
struct add_const_to_root
{
    typedef typename add_const_to_root<T>::type *  type;
};

/// Add const on all levels
/** In case of pointers, constify the address and recursively the type it is pointing to. 
    Other types are constified like in boost::add_const.
    \sa add_const_to_data, add_const_to_all  **/
template <typename T>
struct add_const_to_all
{
    typedef T const  type;
};

template <typename T*>
struct add_const_to_all
{
    typedef typename add_const_to_all<T>::type * const type;
};

}} // namespace mtl::traits

#endif // MTL_TRAITS_ADD_CONST_INCLUDE


