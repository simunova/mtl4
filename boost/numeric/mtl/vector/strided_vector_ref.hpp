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

#ifndef MTL_STRIDED_VECTOR_REF_INCLUDE
#define MTL_STRIDED_VECTOR_REF_INCLUDE


#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <boost/static_assert.hpp>
#include <boost/utility/enable_if.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/common_include.hpp>
#include <boost/numeric/mtl/vector/all_vec_expr.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/vector/crtp_base_vector.hpp>
#include <boost/numeric/mtl/utility/dense_el_cursor.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>
#include <boost/numeric/mtl/utility/strided_dense_el_iterator.hpp>
#include <boost/numeric/mtl/utility/strided_dense_el_cursor.hpp>
#include <boost/numeric/mtl/operation/is_negative.hpp>


namespace mtl { namespace vec {


/// Class for referring vectors stored in strides, e.g. columns in a row-major matrix
template <class Value, typename Parameters = parameters<> >
class strided_vector_ref
  : public vec_expr<strided_vector_ref<Value, Parameters> >,
    public crtp_base_vector< strided_vector_ref<Value, Parameters>, Value, std::size_t >
{
    typedef strided_vector_ref                                                       self;
    typedef crtp_base_vector< self, Value, std::size_t >                             crtp_base;
    typedef crtp_vector_assign< self, Value, std::size_t >                           assign_base;
    typedef vec_expr<strided_vector_ref<Value, Parameters> >                         expr_base;
  public:
    typedef typename boost::remove_const<Value>::type                                value_type ; 
    typedef std::size_t                                                              size_type ;
    typedef Value&                                                                   reference ;
    typedef Value const&                                                             const_reference ;
    typedef Value*                                                                   pointer ;
    typedef Value const*                                                             const_pointer ;
    typedef const_pointer                                                            key_type;
    typedef mtl::strided_dense_el_cursor<Value>                                      cursor_type;
    typedef mtl::strided_dense_el_const_iterator<Value>                              const_iterator;
    typedef mtl::strided_dense_el_iterator<Value>                                    iterator;
    typedef typename Parameters::orientation                                         orientation;
    
    void check_index( size_type MTL_DEBUG_ARG(i) ) const
    {
	MTL_DEBUG_THROW_IF( is_negative(i) || i >= size(*this), index_out_of_range());
    }

    void check_dim( size_type MTL_DEBUG_ARG(s) ) const
    {
	MTL_DEBUG_THROW_IF( size(*this) != 0 && size(*this) != s, incompatible_size());
    }

    template <class E>
    void check_consistent_shape( vec_expr<E> const& ) const
    {
	MTL_DEBUG_THROW_IF((!boost::is_same<
			        typename ashape::ashape<self>::type
			      , typename ashape::ashape<E>::type
			    >::value),
			   incompatible_shape());
    }

  private:
    /// Make default constructor invisible
    strided_vector_ref();

  public:

    /// Constructor take address, length and stride
    strided_vector_ref( size_type length, pointer start_address, size_type stride= 1)
      : data(start_address), my_size(length), my_stride(stride), cloned(false) {}

    /// Clone constructor
    strided_vector_ref(const strided_vector_ref& other, clone_ctor)
      : data(new Value[other.my_size]), my_size(other.my_size), my_stride(1), cloned(true)
    {
	*this= other;
    }

    // Default copy constructor refers to same vector which is okay

    ~strided_vector_ref() { if (cloned && data) delete[] data; }

    //  friend size_type inline size(const self& v) { return v.my_size; } // impedes explicit namespace qualification

    template <typename V2, typename P2>
    friend std::size_t size(const strided_vector_ref<V2, P2>& v);

    size_type stride() const { return my_stride ; }

    reference operator()( size_type i ) { check_index(i); return data[i * my_stride]; }
    const_reference operator()( size_type i ) const { check_index(i); return data[i * my_stride]; }

    reference operator[]( size_type i ) { return (*this)( i ) ; }
    const_reference operator[]( size_type i ) const { return (*this)( i ) ;  }

    self operator[]( irange r ) { return sub_vector(*this, r.start(), r.finish()); }
    const self  operator[]( irange r ) const { return sub_vector(*this, r.start(), r.finish());  }
    
    void delay_assign() const {}

    const_iterator begin() const { return const_iterator(data, my_stride); }
    const_iterator end() const { return const_iterator(data + my_size * my_stride, my_stride); }

    iterator begin() { return iterator(data, my_stride); }
    iterator end() { return iterator(data + my_size * my_stride, my_stride); }

    /// Address of first data entry; to be used with care.
    pointer address_data() { return data; }
    const_pointer address_data() const { return data; }

    // from pointer to index
    size_type offset(const_pointer p) const 
    { 
	size_type o= p - data, i= o / my_stride;
	MTL_DEBUG_THROW_IF(o % my_stride, logic_error("Address not consistent with stride."));
	check_index(i);
	return i;
    }
    
    friend size_type inline num_rows(const self& v) { return mtl::traits::is_row_major<self>::value ? 1 : size(v); }
    friend size_type inline num_cols(const self& v) { return mtl::traits::is_row_major<self>::value ? size(v) : 1; }
    
    vec_vec_asgn_expr<self, self> operator=( self const& e ) 
    {
	return vec_vec_asgn_expr<self, self>( *this, e );
    }

    // self& operator=(self src) Cannot move!

    using assign_base::operator=;

    template <typename Value2> friend void fill(self& vector, const Value2& value)
    {
	std::fill(vector.begin(), vector.end(), value);
    }

    /// Swapping not efficient since elements have to be swapped for not owning the data
    friend void swap(self& vector1, self& vector2)
    {
	vector1.check_dim(vector2.size()); // size(vector2) doesn't compiled with ICC 10.1
	for (size_type i= 0; i < vector1.size(); ++i)
	    swap(vector1[i], vector2[i]);
    }

    void change_dim(size_type MTL_DEBUG_ARG(n)) { MTL_DEBUG_THROW_IF(my_size != n, incompatible_size()); }
    void crop() {} // Only dummy here

  private:
    pointer     data;
    size_type   my_size, my_stride;
    bool        cloned;
} ; // strided_vector_ref


template <typename Value, typename Parameters>
inline std::size_t size(const strided_vector_ref<Value, Parameters>& v)
{
    return v.my_size;
}


template <typename Value, typename Parameters>
strided_vector_ref<Value, Parameters>
inline sub_vector(strided_vector_ref<Value, Parameters>& v, 
		  typename strided_vector_ref<Value, Parameters>::size_type start,
		  typename strided_vector_ref<Value, Parameters>::size_type finish)
{
    using std::min;
    typedef strided_vector_ref<Value, Parameters>    Vector;
    typedef typename Vector::size_type               size_type;

    MTL_DEBUG_THROW_IF( is_negative(start) || is_negative(finish), index_out_of_range());
    finish= min(finish, size(v));
    start= min(start, finish); // implies min(start, size(v))
    return Vector(start <= finish ? finish - start : size_type(0), &v[start], v.stride());
}

template <typename Value, typename Parameters>
const strided_vector_ref<Value, Parameters>
inline sub_vector(const strided_vector_ref<Value, Parameters>& v, 
		  typename strided_vector_ref<Value, Parameters>::size_type start,
		  typename strided_vector_ref<Value, Parameters>::size_type finish)
{
    typedef strided_vector_ref<Value, Parameters>    Vector;
    return sub_vector(const_cast<Vector&>(v), start, finish);
}


}} // namespace mtl::vector

namespace mtl {

    // Enable cloning of strided_vector_ref
    template <typename Value, typename Parameters>
    struct is_clonable< vec::strided_vector_ref<Value, Parameters> > : boost::mpl::true_ {};
        
} // namespace mtl

namespace mtl { namespace traits {


// ================
// Range generators
// For cursors
// ================

    template <typename Value, class Parameters>
    struct range_generator<tag::all, vec::strided_vector_ref<Value, Parameters> >
      : public detail::strided_element_range_generator<
	  vec::strided_vector_ref<Value, Parameters>,
	  const vec::strided_vector_ref<Value, Parameters>,
	  mtl::strided_dense_el_cursor<Value>
	> {};

    template <typename Value, class Parameters>
    struct range_generator<tag::nz, vec::strided_vector_ref<Value, Parameters> >
      : public range_generator<tag::all, vec::strided_vector_ref<Value, Parameters> > {};

    template <typename Value, class Parameters>
    struct range_generator<tag::iter::all, vec::strided_vector_ref<Value, Parameters> >
      : public detail::strided_element_range_generator<
	  vec::strided_vector_ref<Value, Parameters>,
	  vec::strided_vector_ref<Value, Parameters>,
	  mtl::strided_dense_el_iterator<Value>
	> {};

    template <typename Value, class Parameters>
    struct range_generator<tag::iter::nz, vec::strided_vector_ref<Value, Parameters> >
      : public range_generator<tag::iter::all, vec::strided_vector_ref<Value, Parameters> > {};

    template <typename Value, class Parameters>
    struct range_generator<tag::const_iter::all, vec::strided_vector_ref<Value, Parameters> >
      : public detail::strided_element_range_generator<
	  vec::strided_vector_ref<Value, Parameters>,
	  const vec::strided_vector_ref<Value, Parameters>,
	  mtl::strided_dense_el_const_iterator<Value>
	> {};

    template <typename Value, class Parameters>
    struct range_generator<tag::const_iter::nz, vec::strided_vector_ref<Value, Parameters> >
	: public range_generator<tag::const_iter::all, vec::strided_vector_ref<Value, Parameters> >
    {};

	
}} // namespace mtl::traits


#endif // MTL_STRIDED_VECTOR_REF_INCLUDE

