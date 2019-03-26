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

#ifndef MTL_MOD_N_INCLUDE
#define MTL_MOD_N_INCLUDE

#include <iostream>
#include <boost/operators.hpp>
#include <cassert>

#include <boost/config/concept_macros.hpp> 
#ifdef __GXX_CONCEPTS__
#  include <bits/concepts.h>
#endif

#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/linear_algebra/is_invertible.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>
#include <boost/numeric/linear_algebra/operators.hpp>
#include <boost/numeric/linear_algebra/concepts.hpp>
#include <boost/numeric/meta_math/is_prime.hpp>

namespace mtl {

template<typename T, T N>
//  where std::Integral<T>
class mod_n_t 
  : boost::totally_ordered< mod_n_t<T, N> >,
    boost::arithmetic< mod_n_t<T, N> >
{
    // or BOOST_STATIC_ASSERT((IS_INTEGRAL))

    T                 value;
 public:
    typedef T         value_type;
    typedef mod_n_t   self;

    static T const modulo= N;

    mod_n_t() : value(0) {}

    explicit mod_n_t(T const& v)
    {
	value= v >= 0 ? v%modulo : modulo - -v%modulo; 
    }

    // modulo of negative numbers can be bizarre
    // better use constructor for T
    explicit mod_n_t(int v) 
    {
	value= v >= 0 ? v%modulo : modulo - -v%modulo; 
    }

    // copy constructor
    mod_n_t(const mod_n_t<T, N>& m): value(m.get()) {}
  
    // assignment
    mod_n_t<T, N>& operator= (const mod_n_t<T, N>& m) 
    {
	value= m.value; 
	return *this; 
    }

    mod_n_t<T, N>& operator= (T const& v)
    {
	value= v >= 0 ? v%modulo : modulo - -v%modulo; 
	return *this; 
    }
    
    // conversion from other moduli must be called explicitly
    template<T OtherN>
    mod_n_t<T, N>& convert(const mod_n_t<T, OtherN>& m) 
    {
	value= m.value >= 0 ? m.value%modulo : modulo - -m.value%modulo; 
	return *this; 
    }

    T get() const 
    {
	return value; 
    }

    bool operator==(self const& y) const
    {
	check(*this); check(y);
	return this->value == y.value;
    }

    bool operator<(self const& y) const
    {
	check(*this); check(y);
	return this->value < y.value;
    }

    self& operator+= (self const& y)
    {
	check(*this); check(y);
	this->value += y.value;
	this->value %= modulo;
	return *this;
    }

    self& operator-= (self const& y)
    {
	check(*this); check(y);
	// add n to avoid negative numbers esp. if T is unsigned
	this->value += modulo;
	this->value -= y.value;
	this->value %= modulo;
	return *this;
    }

    self& operator*= (self const& y)
    {
	check(*this); check(y);
	this->value *= y.value;
	this->value %= modulo;
	return *this;
    }

    self& operator/= (self const& y);
    
};

template<typename T, T N>
inline void check(const mod_n_t<T, N>& x)
{
    assert(x.get() >= 0 && x.get() < N);
}

template<typename T, T N>
inline std::ostream& operator<< (std::ostream& stream, const mod_n_t<T, N>& a) 
{
    check(a);
    return stream << a.get(); 
}


// Extended Euclidian algorithm in vector notation
    // uu = (u1, u2, u3) := (1, 0, u)
    // vv = (v1, v2, v3) := (0, 1, v)
    // while (u3 % v3 != 0) {
    //   q= u3 / v3
    //   rr= uu - q * vv
    //   uu= vv
    //   vv= rr }
    // 
    // with u = N and v = y
    // --> v2 * v mod u == gcd(u, v)
    // --> v2 * y mod N == 1
    // --> x * v2 == x / y
    // v1, u1, and r1 not used
template<typename T, T N>
inline mod_n_t<T, N>& mod_n_t<T, N>::operator/= (const mod_n_t<T, N>& y) 
{
    check(*this); check(y);
    if (y.get() == 0) throw "Division by 0";

    // Goes wrong with unsigned b/c some values will be negative (even if the result isn't)
    // Something like remove_sign<T>::type would be cute
    int u= N, v= y.get(), /* u1= 1, */  u2= 0, /* v1= 0, */  v2= 1, q, r, /* r1, */  r2;

    while (u % v != 0) {
	q= u / v;

	r= u % v; /* r1= u1 - q * v1; */ r2= u2 - q * v2;
	u= v; /* u1= v1; */ u2= v2;
	v= r; /* v1= r1; */ v2= r2;
    }

    return *this *= mod_n_t<T, N>(v2); 
}

inline int gcd(int u, int v)
{
    int r;
    while ((r= u % v) != 0) {
	u= v; v= r;
    }
    return v;
}

} // namespace mtl

namespace math {

    using mtl::mod_n_t;

    template<typename T, T N>
    struct identity_t< add< mod_n_t<T, N> >, mod_n_t<T, N> >
    {
	typedef mod_n_t<T, N> mod_t;
	mod_t operator() (add<mod_t> const&, mod_t const& v) const
	{
	    return mod_t(0);
	}
    };


    // Reverse definition, a little more efficient if / uses inverse
    template<typename T, T N>
    struct inverse_t< add< mod_n_t<T, N> >, mod_n_t<T, N> >
    {
	typedef mod_n_t<T, N> mod_t;
	mod_t operator() (add<mod_t> const& op, mod_t const& v) const
	{
	    return identity(op, v) - v;
	}
    };
    

    template<typename T, T N>
    struct is_invertible_t< add< mod_n_t<T, N> >, mod_n_t<T, N> >
    {
	typedef mod_n_t<T, N> mod_t;
	bool operator() (add<mod_t> const&, mod_t const& v) const
	{ return true; }
    };
    

    template<typename T, T N>
    struct identity_t< mult< mod_n_t<T, N> >, mod_n_t<T, N> >
    {
	typedef mod_n_t<T, N> mod_t;
	mod_t operator() (mult<mod_t> const&, mod_t const& v) const
	{
	    return mod_t(1);
	}
    };


    // Reverse definition, a little more efficient if / uses inverse
    template<typename T, T N>
    struct inverse_t< mult< mod_n_t<T, N> >, mod_n_t<T, N> >
    {
	typedef mod_n_t<T, N> mod_t;
	mod_t operator() (mult<mod_t> const&, mod_t const& v) const
	{
	    return mod_t(1) / v;
	}
    };
    

    template<typename T, T N>
    struct is_invertible_t< mult< mod_n_t<T, N> >, mod_n_t<T, N> >
    {
	typedef mod_n_t<T, N> mod_t;
	bool operator() (mult<mod_t> const&, mod_t const& v) const
	{
#           ifdef MTL_TRACE_MOD_N_INVERTIBILITY_DISPATCHING 
                std::cout << "[slow mod n inversion test] ";
#           endif
	    T value = v.get();
	    return value != 0 && mtl::gcd(N, value) == 1;
	}
    };
    

# ifdef __GXX_CONCEPTS__

    // With Concept we can provide a faster invertibility test for prime numbers:
    //  only 0 is not invertible and gcd doesn't need to be called

    template<typename T, T N>
        where meta_math::Prime<N>
    struct is_invertible_t< mult< mod_n_t<T, N> >, mod_n_t<T, N> >
    {
	typedef mod_n_t<T, N> mod_t;
	bool operator() (mult<mod_t> const&, mod_t const& v) const
	{
#           ifdef MTL_TRACE_MOD_N_INVERTIBILITY_DISPATCHING 
                std::cout << "[fast mod n inversion test] ";
#           endif
 	    return v.get() != 0;
	}
    };




// Concept mapping
// All modulo sets are commutative rings with identity
// but only if N is prime it is also a field
// Due to some mapping nesting trouble we define normally derived maps


template <typename T, T N>
concept_map CommutativeRingWithIdentity< mod_n_t<T, N> > 
{
    // Why do we need the typedefs???
    
    typedef mod_n_t<T, N>& plus_assign_result_type;
    typedef mod_n_t<T, N>  addition_result_type;
    typedef mod_n_t<T, N>  unary_result_type;
    typedef mod_n_t<T, N>& minus_assign_result_type;
    typedef mod_n_t<T, N>  subtraction_result_type;

    typedef mod_n_t<T, N>& mult_assign_result_type;
    typedef mod_n_t<T, N>  mult_result_type;
    typedef mod_n_t<T, N>& divide_assign_result_type;
    typedef mod_n_t<T, N>  division_result_type;

    typedef mod_n_t<T, N>  inverse_result_type;
    typedef mod_n_t<T, N>  identity_result_type;
    typedef bool           is_invertible_result_type;
}

template <typename T, T N>
concept_map MultiplicativePartiallyInvertibleMonoid< mod_n_t<T, N> >
{
    // Why do we need the typedefs???

    typedef mod_n_t<T, N>& mult_assign_result_type;
    typedef mod_n_t<T, N>  mult_result_type;
    typedef mod_n_t<T, N>& divide_assign_result_type;
    typedef mod_n_t<T, N>  division_result_type;

    typedef mod_n_t<T, N>  inverse_result_type;
    typedef mod_n_t<T, N>  identity_result_type;
    typedef bool           is_invertible_result_type;
}


template <typename T, T N>
    where meta_math::Prime<N>
concept_map Field< mod_n_t<T, N> >
{
    // Why do we need the typedefs???

    typedef mod_n_t<T, N>& plus_assign_result_type;
    typedef mod_n_t<T, N>  addition_result_type;
    typedef mod_n_t<T, N>  unary_result_type;
    typedef mod_n_t<T, N>& minus_assign_result_type;
    typedef mod_n_t<T, N>  subtraction_result_type;

    typedef mod_n_t<T, N>& mult_assign_result_type;
    typedef mod_n_t<T, N>  mult_result_type;
    typedef mod_n_t<T, N>& divide_assign_result_type;
    typedef mod_n_t<T, N>  division_result_type;

    typedef mod_n_t<T, N>  inverse_result_type;
    typedef mod_n_t<T, N>  identity_result_type;
    typedef bool           is_invertible_result_type;
}

# endif // __GXX_CONCEPTS__

} // namespace math


#endif // MTL_MOD_N_INCLUDE
