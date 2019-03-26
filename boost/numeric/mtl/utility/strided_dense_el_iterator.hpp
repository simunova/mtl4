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

#ifndef MTL_STRIDED_DENSE_EL_ITERATOR_INCLUDE
#define MTL_STRIDED_DENSE_EL_ITERATOR_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/detail/strided_base_cursor.hpp>

namespace mtl {


/// Iterator going in strides over element of matrix, matrix row/column, or vector
/** - Strided iterator *operator returns (const) reference to Value instead of key
    - row(i) and col(i) don't work 
**/
template <typename Value>
struct strided_dense_el_const_iterator
    : public detail::strided_base_cursor<const Value*> 
{
    typedef const Value*                              key_type;
    typedef detail::strided_base_cursor<key_type>     super;
    typedef strided_dense_el_const_iterator           self;

    strided_dense_el_const_iterator(key_type me, size_t stride) : super(me, stride) {}

    template <typename Parameters>
    strided_dense_el_const_iterator(mtl::mat::dense2D<Value, Parameters> const& ma, size_t r, size_t c, size_t stride)
	: super(ma.elements() + ma.indexer(ma, r, c), stride)
    {}

    self operator+(int x) const
    {
	return super::operator+(x);
    }

    const Value& operator*() const
    {
	return *(this->key);
    }
};

/// Iterator going in strides over element of matrix, matrix row/column, or vector
/** - Strided iterator *operator returns (const) reference to Value instead of key
    - row(i) and col(i) don't work 
**/
template <typename Value>
struct strided_dense_el_iterator
    : public detail::strided_base_cursor<Value*> 
{
    typedef Value*                                    key_type;
    typedef detail::strided_base_cursor<key_type>     super;
    typedef strided_dense_el_iterator                 self;

    strided_dense_el_iterator(key_type me, size_t stride) : super(me, stride) {}

    template <typename Parameters>
    strided_dense_el_iterator(mtl::mat::dense2D<Value, Parameters>& ma, size_t r, size_t c, size_t stride)
	: super(ma.elements() + ma.indexer(ma, r, c), stride)
    {}

    self operator+(int x) const
    {
	self tmp(*this);
	tmp+= x;
	return tmp;
    }

    Value& operator*() const
    {
	return *(this->key);
    }
};


} // namespace mtl

#endif // MTL_STRIDED_DENSE_EL_ITERATOR_INCLUDE
