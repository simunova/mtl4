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

#ifndef MTL_STRIDED_BASE_CURSOR_INCLUDE
#define MTL_STRIDED_BASE_CURSOR_INCLUDE

#include <boost/numeric/mtl/detail/base_cursor.hpp>

namespace mtl { namespace detail {

template <class Key> struct strided_base_cursor 
 : base_cursor<Key>
{
    typedef Key                  key_type;
    typedef base_cursor<Key>     base;
    typedef strided_base_cursor  self;

    strided_base_cursor () {} 
    strided_base_cursor (key_type key, std::size_t stride) 
	: base(key), stride(stride) 
    {}

    self& operator++ () 
    { 
	this->key+= stride; return *this; 
    }
    self operator++ (int) 
    { 
	self tmp = *this; 
	this->key+= stride; 
	return tmp; 
    }
    self& operator-- () 
    { 
	this->key-= stride; 
	return *this; 
    }
    self operator-- (int) 
    { 
	self tmp = *this; 
	this->key-= stride; 
	return tmp; 
    }
    self& operator+=(int n) 
    { 
	this->key += stride * n; 
	return *this; 
    }
    self operator+(int n) const
    {
	self tmp(*this);
	tmp+= n;
	return tmp;
    }
    self& operator-=(int n) 
    { 
	this->key -= stride * n; 
	return *this; 
    }

    int operator-(const self& cc) const 
    {
	return (this->key - cc.key) / stride;
    }

    std::size_t stride;
};

}} // namespace mtl::detail

#endif // MTL_STRIDED_BASE_CURSOR_INCLUDE
