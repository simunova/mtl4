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

#ifndef MTL_MULTI_TMP_INCLUDE
#define MTL_MULTI_TMP_INCLUDE

namespace mtl {

/// Helper class to define a set of temporaries
template <unsigned Size, typename Value>
struct multi_tmp 
{
    typedef multi_tmp<Size-1, Value> sub_type;
	
    multi_tmp() {}
    multi_tmp(const Value& v) : value(v), sub(v) {}
	
    Value sum() { return value + sub.sum(); }
	
    Value    value;
    sub_type sub;
};

template <typename Value>
struct multi_tmp<0, Value> 
{
    multi_tmp() {}
    multi_tmp(const Value&) {}
    Value sum() { return 0; }
};

/// Helper class to define a set of constants and initialize it from an array
template <unsigned Index, unsigned Size, typename Value>
struct multi_constant_from_array
{
    typedef multi_constant_from_array<Index+1, Size, Value> sub_type;

    template <typename Array, typename IndexType>
    multi_constant_from_array(const Array& x, IndexType i) : value(x[i+Index]), sub(x, i) {}

    const Value value;
    sub_type    sub;
};


template <unsigned Size, typename Value>
struct multi_constant_from_array<Size, Size, Value>
{
    template <typename Array, typename IndexType>
    multi_constant_from_array(const Array&, IndexType) {}
};



} // namespace mtl

#endif // MTL_MULTI_TMP_INCLUDE
