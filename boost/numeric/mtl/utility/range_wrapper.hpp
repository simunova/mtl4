// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_TRAITS_RANGE_WRAPPER_INCLUDE
#define MTL_TRAITS_RANGE_WRAPPER_INCLUDE

#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>

namespace mtl { namespace traits {

/// Wrapper for range_generator
/** Instead of passing the collection or cursor to the begin() and end() function,
    it is passed to the constructor.
    As a consequence, it can be used with ranged for from C++11. **/
template <typename Tag, typename Collection>
struct range_wrapper
  : range_generator<Tag, Collection>
{
    typedef range_generator<Tag, Collection> gen_type;
    MTL_STATIC_ASSERT(gen_type::level != 0, "Template arguments not supported, probably traversal with unsupported tag combination.");
    typedef typename gen_type::type          type;

    explicit range_wrapper(const Collection& c) : c(c) {} //< Initialize with collection
    
    type begin() const { return gen_type::begin(c); } //< nullary begin 
    type end() const { return gen_type::end(c); } //< nullary end
    
  private:
    const Collection& c;
};

}} // namespace mtl::traits

namespace mtl { 

/// Cursor over rows of a collection
template <typename T>
mtl::traits::range_wrapper<tag::row, T> 
inline rows_of(const T& x)
{
    return mtl::traits::range_wrapper<tag::row, T>(x);
}

/// Cursor over cols of a collection
template <typename T>
mtl::traits::range_wrapper<tag::col, T> 
inline cols_of(const T& x)
{
    return mtl::traits::range_wrapper<tag::col, T>(x);
}

/// Cursor over major index of a collection
template <typename T>
mtl::traits::range_wrapper<tag::major, T> 
inline major_of(const T& x)
{
    return mtl::traits::range_wrapper<tag::major, T>(x);
}

/// Cursor over minor index of a collection
template <typename T>
mtl::traits::range_wrapper<tag::minor, T> 
inline minor_of(const T& x)
{
    return mtl::traits::range_wrapper<tag::minor, T>(x);
}

/// Cursor over non-zero elements of a collection
template <typename T>
mtl::traits::range_wrapper<tag::nz, T> 
inline nz_of(const T& x)
{
    return mtl::traits::range_wrapper<tag::nz, T>(x);
}

/// Cursor over all elements of a collection
template <typename T>
mtl::traits::range_wrapper<tag::all, T> 
inline all_of(const T& x)
{
    return mtl::traits::range_wrapper<tag::all, T>(x);
}

/// Cursor over tagged range of a collection
template <typename Tag, typename T>
mtl::traits::range_wrapper<Tag, T> 
inline range_of(const T& x)
{
    return mtl::traits::range_wrapper<Tag, T>(x);
}


namespace mat {
    using mtl::rows_of;
    using mtl::cols_of;
    using mtl::major_of;
    using mtl::minor_of;
    using mtl::nz_of;
    using mtl::all_of;
    using mtl::range_of;
} // matrix




} // mtl


#endif // MTL_TRAITS_RANGE_WRAPPER_INCLUDE
