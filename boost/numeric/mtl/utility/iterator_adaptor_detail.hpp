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

#ifndef MTL_ITERATOR_ADAPTOR_DETAIL_INCLUDE
#define MTL_ITERATOR_ADAPTOR_DETAIL_INCLUDE

namespace mtl { namespace utilities { namespace detail {


template <typename Adaptor>
struct adaptor_operators
{
    Adaptor& operator++() 
    {
	Adaptor& me = static_cast<Adaptor&>(*this);
	++me.cursor;
	return me;
    }

    Adaptor& operator++(int) 
    {
	Adaptor& me = static_cast<Adaptor&>(*this);
	Adaptor  tmp(me);
	++me.cursor;
	return tmp;
    }
    
    bool operator==(Adaptor const& x) const
    {
	Adaptor const& me = static_cast<Adaptor const&>(*this);

	// Sloppy, nothing tested about property map
	return me.cursor == x.cursor;

	// Compare addresses of property maps
	// Problem: different addresses doesn't imply that the maps are different
	// return &me.map == &x.map && me.cursor == x.cursor;

	// Certainly better they provide comparison
	// Problem: not guaranteed to exist or can be ridiculously expensive
	// return me.map == x.map && me.cursor == x.cursor; 
    }

    bool operator!=(Adaptor const& x) const
    {
	return !operator==(x);
    }
};

template <typename Adaptor>
struct ra_adaptor_operators 
  : public adaptor_operators<Adaptor>
{
    Adaptor operator+(int i) const
    {
	const Adaptor& me = static_cast<const Adaptor&>(*this);
	return Adaptor(me.map, me.cursor + i);
    }

    Adaptor& operator+=(int i)
    {
	Adaptor& me = static_cast<Adaptor&>(*this);
	me.cursor+= i;
	return me;
    }
};




template <typename PropertyMap, typename Cursor, typename ValueType>
struct const_iterator_proxy
{
    const_iterator_proxy(PropertyMap map, Cursor cursor) 
	: map(map), cursor(cursor) {}

    operator ValueType() const
    {
	return map(*cursor);
    }

    PropertyMap            map;
    Cursor                 cursor;
};


template <typename PropertyMap, typename Cursor, typename ValueType>
struct iterator_proxy
{
    typedef iterator_proxy                    self;

    iterator_proxy(PropertyMap map, Cursor cursor) 
	: map(map), cursor(cursor) {}

    operator ValueType() const
    {
	return map(*cursor);
    }

    self& operator=(ValueType const& value)
    {
	map(*cursor, value);
	return *this;
    }

    PropertyMap           map;
    Cursor                 cursor;
};

}}} // namespace mtl::utilities::detail

#endif // MTL_ITERATOR_ADAPTOR_DETAIL_INCLUDE
