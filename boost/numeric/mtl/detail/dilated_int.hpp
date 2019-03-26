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
//
// Written by Jiahu Deng and Peter Gottschling

#ifndef MTL_DILATED_INT_INCLUDE
#define MTL_DILATED_INT_INCLUDE

#include <boost/numeric/mtl/detail/masked_dilation_tables.hpp>
#include <iostream>

namespace mtl { namespace dilated {

template <typename T>
struct even_bits
{
    static T const value = T(-1) / T(3);
};
    
template <typename T>
struct odd_bits
{
    static T const value = ~even_bits<T>::value;
};

// And is mostly used with original mask
template <typename T, T BitMask, bool Normalized>
struct masking
{
    inline void operator() (T& x) const
    {
	x &= BitMask;
    }
};

// Or is mostly used with complementary mask
template <typename T, T BitMask>
struct masking<T, BitMask, false>
{
    inline void operator() (T& x) const
    {
	static T const anti_mask = ~BitMask;
	x |= anti_mask;
    }    
};

template <typename T, T BitMask>
struct last_bit;

template <typename T, T BitMask, bool IsZero>
struct last_bit_helper {
    static T const tmp = BitMask >> 1;
    static T const value = BitMask & 1 ? 1 : last_bit_helper<T, tmp, tmp == 0>::value << 1;
};

template <typename T, T BitMask>
struct last_bit_helper<T, BitMask, true> {
    static T const value = 0;
};

template <typename T, T BitMask>
struct last_bit
{
    static T const value = last_bit_helper<T, BitMask, BitMask == 0>::value;
};


template <typename T, T BitMask, bool Normalized>
struct dilated_int
{
    typedef T                                       value_type;
    typedef dilated_int<T, BitMask, Normalized>     self;
    
    typedef masking<T, BitMask, Normalized>         clean_carry;
    typedef masking<T, BitMask, !Normalized>        init_carry;

    static T const       bit_mask = BitMask,
	                 anti_mask = ~BitMask,    
	                 dilated_zero = Normalized ? 0 : anti_mask,
			 dilated_one = dilated_zero + last_bit<T, bit_mask>::value;
protected:
    // masked_dilation_tables<T, bit_mask>   mask_tables;
    // masked_dilation_tables<T, anti_mask>  anti_tables; probably not needed
      
// will be protected later
public:              
    T i;

    void dilate(T x)
    {
		static const T to_switch_on = Normalized ? 0 : anti_mask;
		// auto aa = mask<bit_mask>(x) | to_switch_on;
		// auto ab = mask<static_cast<T>(bit_mask)>(x) | to_switch_on
		i = mask<bit_mask>(x) | to_switch_on;
    }

public:

    // Default constructor
    dilated_int()
    {
	i = Normalized ? 0 : anti_mask;
    }
    
    // Only works for odd and even bits and 4-byte-int at this point !!!!!!!!!!!!!!!!!!!
    explicit dilated_int(T x)
    {
	dilate(x);
    }

    // Only works for odd and even bits and 4-byte-int at this point !!!!!!!!!!!!!!!!!!!
    T undilate()
    {
	return unmask<bit_mask>(i);
    }

    T dilated_value() const
    {
	return i;
    }

    self& operator= (self const& x)
    {
	i = x.i;
	return *this;
    }

    self& operator= (T x)
    {
	dilate(x);
	return *this;
    }

    self& operator++ ()
    {
	static T const x = Normalized ? bit_mask : T(-1);
	i -= x;
	clean_carry()(i);
	return *this;
    }

    self operator++ (int)
    {
	self tmp(*this);
	++*this;
	return tmp;
    }

    self& operator+= (self const& x)
    {
	init_carry()(i);
	i+= x.i;
	clean_carry()(i);
	return *this;
    }

    self operator+ (self const& x)
    {
	self tmp(*this);
	return tmp += x;
    }

    self& operator-- ()
    { 
	i -= dilated_one;
	clean_carry()(i);
	return *this;
    }

    self operator-- (int)
    {
	self tmp(*this);
	--*this;
	return tmp;
    }

    self& operator-= (self const& x)
    {
	i -= x.i;
	clean_carry()(i);
	return *this;
    }
	
    self operator- (self const& x) const
    {
	self tmp(*this);
	return tmp -= x;
    }
    
    // advance in both directions, special care is needed for negative values
    self& advance(long inc)
    {
	value_type incv(inc >= 0 ? inc : -inc);
	self incd(incv);
	if (inc >= 0)
	    operator+=(incd);
	else
	    operator-=(incd);
	return *this;
    }

    bool operator== (self const& x) const
    {
	return i == x.i;
    }

    bool operator!= (self const& x) const
    {
	return i != x.i;
    }

    bool operator<= (self const& x) const
    {
	return i <= x.i;
    }

    bool operator< (self const& x) const
    {
	return i < x.i;
    }

    bool operator>= (self const& x) const
    {
	return i >= x.i;
    }

    bool operator> (self const& x) const
    {
	return i > x.i;
    }


};

} // namespace mtl::dilated

using dilated::dilated_int;

} // namespace mtl

template <typename T, T BitMask, bool Normalized>
inline std::ostream& operator<< (std::ostream& os, mtl::dilated::dilated_int<T, BitMask, Normalized> d)
{
    os.setf(std::ios_base::hex, std::ios_base::basefield);
    return os << d.i;
}

#endif // MTL_DILATED_INT_INCLUDE
