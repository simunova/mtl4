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

#ifndef MTL_ZIPPED_SORT_INCLUDE
#define MTL_ZIPPED_SORT_INCLUDE

// Designed for pointers so far and not tested for general iterators
// For internal use only

#include <cmath>
#include <utility>
#include <iterator>
#include <boost/numeric/mtl/utility/exception.hpp>

namespace mtl { namespace utility {

template <typename T, typename U> struct zip_ref;
template <typename T, typename U> struct zip_it;
template <typename T, typename U> struct zip_value;

struct less_0
{
    template <typename T, typename U>
    bool operator()(const zip_ref<T, U>& x, const zip_ref<T, U>& y) const
    {
	return x.a[x.p] < y.a[y.p];
    }
    
    template <typename T, typename U>
    bool operator()(const zip_ref<T, U>& x, const zip_value<T, U>& y) const
    {
	return x.a[x.p] < y.x;
    }
    
    template <typename T, typename U>
    bool operator()(const zip_value<T, U>& x, const zip_ref<T, U>& y) const
    {
	return x.x < y.a[y.p];
    }
};

struct greater_0
{
    template <typename T, typename U>
    bool operator()(const zip_ref<T, U>& x, const zip_ref<T, U>& y) const
    {
	return x.a[x.p] > y.a[y.p];
    }
    
    template <typename T, typename U>
    bool operator()(const zip_ref<T, U>& x, const zip_value<T, U>& y) const
    {
	return x.a[x.p] > y.x;
    }
    
    template <typename T, typename U>
    bool operator()(const zip_value<T, U>& x, const zip_ref<T, U>& y) const
    {
	return x.x > y.a[y.p];
    }
};

struct abs_greater_0
{
    template <typename T, typename U>
    bool operator()(const zip_ref<T, U>& x, const zip_ref<T, U>& y) const
    {
	using std::abs;
	return abs(x.a[x.p]) > abs(y.a[y.p]);
    }
    
    template <typename T, typename U>
    bool operator()(const zip_ref<T, U>& x, const zip_value<T, U>& y) const
    {
	using std::abs;
	return abs(x.a[x.p]) > abs(y.x);
    }
    
    template <typename T, typename U>
    bool operator()(const zip_value<T, U>& x, const zip_ref<T, U>& y) const
    {
	using std::abs;
	return abs(x.x) > abs(y.a[y.p]);
    }
};


template <typename T, typename U>
struct zip_it
{
    typedef zip_ref<T, U>        ref_type;
    typedef long                 diff_type;
    // typedef std::difference_type diff_type;

    explicit zip_it(T* a, U* v, std::size_t p) : a(a), v(v), p(diff_type(p)) {}

    ref_type operator*() { return ref_type(a, v, p); }
    zip_it& operator++() { p++; return *this;}
    zip_it operator++(int) { zip_it tmp(a, v, p); p++; return tmp;}
    zip_it& operator--() { p--; return *this;}
    zip_it operator--(int) { zip_it tmp(a, v, p); p--; return tmp;}

    void check(const zip_it& MTL_DEBUG_ARG(other)) const { assert(a == other.a); assert(v == other.v); }

    bool operator==(const zip_it& other) const { check(other); return p == other.p; }
    bool operator!=(const zip_it& other) const { check(other); return p != other.p; }
    bool operator<=(const zip_it& other) const { check(other); return p <= other.p; }
    bool operator<(const zip_it& other) const { check(other); return p < other.p; }
    bool operator>=(const zip_it& other) const { check(other); return p >= other.p; }
    bool operator>(const zip_it& other) const { check(other); return p > other.p; }
    diff_type operator-(const zip_it& other) const { check(other); return p - other.p; }
    template <typename S> zip_it operator+(S i) const { return zip_it(a, v, p+i); } // use template to avoid narrowing warnings
    template <typename S> zip_it& operator+=(S i) { p += diff_type(i); return *this; } // dito
    template <typename S> zip_it operator-(S i) const { return zip_it(a, v, p - diff_type(i)); } // dito
    zip_it& operator=(const zip_it& other) { check(other); p= other.p; return *this; }

    T*            a;
    U*            v;
    diff_type     p;
};


template <typename T, typename U>
struct zip_ref
{
    typedef zip_ref       self;

    zip_ref(T* a, U* v, int p) : a(a), v(v), p(p) {}

    void check(const zip_ref& MTL_DEBUG_ARG(other)) const { assert(a == other.a); assert(v == other.v); }

    bool operator<(const zip_ref& r) const { check(r); return a[p] < r.a[r.p]; }
    zip_ref& operator=(const zip_ref& r) 
    { 
	check(r);
	if (p == r.p) 
	    return *this;
	a[p]= r.a[r.p];	v[p]= r.v[r.p];
	p= r.p; 
	return *this;
    }
#ifdef MTL_WITH_MOVE
    zip_ref& operator=(zip_value<T, U>&& zv) 
    {
	a[p]= zv.x;
	v[p]= zv.y;
	return *this;
    }
#endif

    zip_ref& operator=(const zip_value<T, U>& zv)
    {
	a[p]= zv.x;
	v[p]= zv.y;
	return *this;
    }

    T *a;
    U *v;
    int           p;

};

// const ref is ugly but mutable ref doesn't work with temporaries and copy is not as good as overload
template <typename T, typename U>
inline void swap(const zip_ref<T, U>& x, const zip_ref<T, U>& y)
{
    using std::swap;

    swap(x.a[x.p], y.a[y.p]);
    swap(x.v[x.p], y.v[y.p]);
}    

template <typename T, typename U>
struct zip_value
{
    zip_value(const zip_ref<T, U>& r) : x(r.a[r.p]), y(r.v[r.p]) {}

    T x;
    U y;
};


}} // namespace mtl::utility

namespace std {
    template <typename T, typename U>
    struct iterator_traits<mtl::utility::zip_it<T, U> >
    {
	typedef mtl::utility::zip_ref<T, U>    ref_type;
	typedef mtl::utility::zip_value<T, U>  value_type;
	typedef ref_type&             reference;
	typedef ref_type*             pointer;
	typedef int                   difference_type;
	typedef random_access_iterator_tag iterator_category;
    };
}

// usage:
// sort(zip_it<T, U>(a, v, 0), zip_it<T, U>(a, v, S), less_0());
// where a and v are pointers or arrays and S the size of both arrays

#endif // MTL_ZIPPED_SORT_INCLUDE
