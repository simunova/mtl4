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

#ifndef MTL_STRIDED_DENSE_EL_CURSOR_INCLUDE
#define MTL_STRIDED_DENSE_EL_CURSOR_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/detail/strided_base_cursor.hpp>

namespace mtl {

/// Cursor going in strides over element of matrix, matrix row/column, or vector
template <typename Value> 
struct strided_dense_el_cursor : public detail::strided_base_cursor<const Value*> 
{
    typedef Value                                        value_type;
    typedef const value_type*                            const_pointer_type;
    typedef detail::strided_base_cursor<const Value*>    super;
    typedef strided_dense_el_cursor                      self;

    //  strided_dense_el_cursor () {} 
    strided_dense_el_cursor (const_pointer_type me, size_t stride) : super(me, stride) {}

    template <typename Parameters>
    strided_dense_el_cursor(mtl::mat::dense2D<Value, Parameters> const& ma, size_t r, size_t c, size_t stride)
	: super(ma.elements() + ma.indexer(ma, r, c), stride)
    {}

    // Why do we need this?
    strided_dense_el_cursor(super const& x) : super(x) {}

    self operator+(int x) const
    {
	return super::operator+(x);
    }
};

} // namespace mtl

#endif // MTL_STRIDED_DENSE_EL_CURSOR_INCLUDE
