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

#ifndef MTL_MATRIX_MULTI_VECTOR_INCLUDE
#define MTL_MATRIX_MULTI_VECTOR_INCLUDE

#include <cassert>
#include <boost/utility/enable_if.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/matrix/crtp_base_matrix.hpp>
#include <boost/numeric/mtl/matrix/mat_expr.hpp>
#include <boost/numeric/mtl/matrix/multi_vector_range.hpp>
#include <boost/numeric/mtl/matrix/make_fast_multi_vector_expr.hpp>
#include <boost/numeric/mtl/utility/is_what.hpp>
#include <boost/numeric/mtl/utility/is_composable_vector.hpp>
#include <boost/numeric/mtl/utility/is_multi_vector_expr.hpp>
#include <boost/numeric/mtl/utility/fast_multi_vector_expr.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>


namespace mtl { namespace mat {

    
// Might need to be defined later
struct multi_vector_key {};

/// Matrix constituting of set of column vectors (under development)
template <typename Vector>
class multi_vector
  : public base_matrix<typename mtl::Collection<Vector>::value_type, parameters<> >,
    public crtp_base_matrix< multi_vector<Vector>, typename Collection<Vector>::value_type, 
			     typename Collection<Vector>::size_type>,
    public mat_expr< multi_vector<Vector> >    
{
    typedef base_matrix<typename Collection<Vector>::value_type, parameters<> >           super;

    // Vector must by column vector
    MTL_STATIC_ASSERT((boost::is_same<typename OrientedCollection<Vector>::orientation,
			              tag::col_major>::value),
		      "Vector must be a column vector.");
  public:
    typedef multi_vector                             self;
    // typedef mtl::mat::parameters<>                parameters;
    typedef Vector                                   vector_type;
    typedef tag::col_major                           orientation;
    typedef typename Collection<Vector>::value_type  value_type;
    typedef typename Collection<Vector>::size_type   size_type;
    typedef const value_type&                        const_reference;
    typedef value_type&                              reference;
    typedef multi_vector_key                         key_type;
    typedef crtp_matrix_assign< self, value_type, size_type >    assign_base;

    multi_vector() : super(non_fixed::dimensions(0, 0)), master(0) {}

  private:
    void setup_data(size_type num_rows, size_type num_cols, boost::mpl::false_)
    {
	for (size_type i= 0; i < num_cols; ++i)
	    data[i]= Vector(num_rows);
	master= 0;
    }
    
    void setup_data(size_type num_rows, size_type num_cols, boost::mpl::true_)
    {
	master= new Vector(num_rows * num_cols);
	for (size_type i= 0; i < num_cols; ++i) {
	    Vector tmp((*master)[irange(i * num_rows, (i+1) * num_rows)]);
	    swap(data[i], tmp);
	}
    }

  public:
    /// Constructor by number of rows and columns
    multi_vector(size_type num_rows, size_type num_cols)
      : super(non_fixed::dimensions(num_rows, num_cols)), data(num_cols)
    {
	setup_data(num_rows, num_cols, mtl::traits::is_composable_vector<Vector>());
	this->my_nnz= num_rows * num_cols;
    }

    /// Constructor column vector and number of columns (for easier initialization)
    multi_vector(const Vector& v, size_type num_cols)
      : super(non_fixed::dimensions(size(v), num_cols)), data(num_cols)
    {
	using mtl::num_rows;
	setup_data(num_rows(v), num_cols, mtl::traits::is_composable_vector<Vector>());
	for (size_type i= 0; i < num_cols; ++i)
	    data[i]= v;
	this->my_nnz= num_cols * size(v);
    }

    ~multi_vector() { delete master; }

  private:

    void change_dim(size_type r, size_type c, boost::mpl::false_)
    {
	for (size_type i= 0; i < c; i++)
	    data[i].change_dim(r);
    }

    void change_dim(size_type r, size_type c, boost::mpl::true_)
    {
	delete master;
	setup_data(r, c, boost::mpl::true_());
    }

  public:
    /// Change dimension, can keep old data
    void change_dim(size_type r, size_type c)
    {
	if (r == super::num_rows() && c == super::num_cols()) return; 
	super::change_dim(r, c);
	data.change_dim(c);
	change_dim(r, c, mtl::traits::is_composable_vector<Vector>());
    }

    void self_assignment(const self& src, boost::mpl::false_)
    {
	for (size_type i= 0; i < super::num_cols(); i++)
	    data[i]= src.data[i];
    }

    void self_assignment(const self& src, boost::mpl::true_)
    {	*master= *src.master;    }

    /// Copy constructor 
    /** Explicitly needed now. **/
    self& operator=(const self& src)
    {
	assign_base::checked_change_dim(src.num_rows(), src.num_cols());
	self_assignment(src, mtl::traits::is_composable_vector<Vector>());
	return *this;
    }

    // Todo: multi_vector with other matrix expressions
    /// Assign multi_vector and expressions thereof, general matrices currently not allowed 
    template <typename Src>
    typename boost::enable_if_c<mtl::traits::is_multi_vector_expr<Src>::value
				&& !mtl::traits::is_fast_multi_vector_expr<Src>::value, self&>::type
    operator=(const Src& src)
    {
	MTL_THROW_IF((mtl::mat::num_rows(src) != super::num_rows() 
		      || mtl::mat::num_cols(src) != super::num_cols()), incompatible_size());
	for (std::size_t i= 0, n= super::num_cols(); i < n; ++i)
	    vector(i)= src.vector(i);
	return *this;
    }

    /// Assign multi_vector and expressions thereof that are stored internally by a single vector
    template <typename Src>
    typename boost::enable_if<mtl::traits::is_fast_multi_vector_expr<Src>, self&>::type
    operator=(const Src& src)
    {
	assert(master);
	*master= make_fast_multi_vector_expr(src);
	return *this;
    }

    template <typename Src>
    typename boost::enable_if_c<mtl::traits::is_matrix<Src>::value 
				&& !mtl::traits::is_multi_vector_expr<Src>::value, self&>::type
    operator=(const Src& src)
    {
	assign_base::operator=(src);
	return *this;
    }
    
    /// Assign scalar
    template <typename Src>
    typename boost::enable_if<mtl::traits::is_scalar<Src>, self&>::type
    operator=(const Src& src)
    {
	assign_base::operator=(src);
	return *this;
    }

    const_reference operator() (size_type i, size_type j) const { return data[j][i]; }
    reference operator() (size_type i, size_type j) { return data[j][i]; }

    Vector& vector(size_type i) { return data[i]; }
    const Vector& vector(size_type i) const { return data[i]; }

    Vector& at(size_type i) { return data[i]; }
    const Vector& at(size_type i) const { return data[i]; }

    multi_vector_range<Vector> vector(irange const& r) const { return multi_vector_range<Vector>(*this, r); }

    inline friend const Vector& make_fast_multi_vector_expr(const self& mv) { assert(mv.master); return *mv.master; }
    inline friend Vector& make_fast_multi_vector_expr(self& mv) { assert(mv.master); return *mv.master; }

  protected:  
    mtl::dense_vector<Vector, mtl::vec::parameters<> >          data;
    Vector*                                                                master;
};

/// Number of rows
template< typename Vector >
typename Collection< Vector >::size_type num_cols(const multi_vector< Vector >& A) { return A.num_cols(); }

/// Number of columns
template< typename Vector >
typename Collection< Vector >::size_type num_rows(const multi_vector< Vector >& A) { return A.num_rows(); }

/// Size as defined by number of rows times columns
template< typename Vector >
typename Collection< Vector >::size_type size(const multi_vector< Vector >& A) { return num_rows(A) * num_cols(A); }
}} // namespace mtl::matrix

namespace mtl {
	using mat::multi_vector;
}

#endif // MTL_MATRIX_MULTI_VECTOR_INCLUDE
