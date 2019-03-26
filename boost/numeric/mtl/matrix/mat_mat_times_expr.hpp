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

#ifndef MTL_MAT_MAT_TIMES_EXPR_INCLUDE
#define MTL_MAT_MAT_TIMES_EXPR_INCLUDE

#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/type_traits/is_base_of.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/matrix/mat_mat_op_expr.hpp>
#include <boost/numeric/mtl/operation/sfunctor.hpp>
#include <boost/numeric/mtl/operation/compute_factors.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/concept/std_concept.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>


namespace mtl { namespace mat {

template <typename E1, typename E2>
struct mat_mat_times_expr 
  : public mat_mat_op_expr< E1, E2, mtl::sfunctor::times<typename Collection<E1>::value_type, typename Collection<E2>::value_type> >,
    public mat_expr< mat_mat_times_expr<E1, E2> >
{
    typedef mat_mat_op_expr< E1, E2, mtl::sfunctor::times<typename Collection<E1>::value_type, typename Collection<E2>::value_type> > op_base;
    typedef mat_expr< mat_mat_times_expr<E1, E2> >                                                       crtp_base;
    typedef mat_mat_times_expr                   self;
    typedef E1                                   first_argument_type ;
    typedef E2                                   second_argument_type ;
    typedef typename E1::orientation             orientation;
    typedef mtl::non_fixed::dimensions           dim_type;
    typedef typename E1::key_type                key_type;


    typedef typename Collection<E1>::value_type  first_value_type;
    typedef typename Collection<E2>::value_type  second_value_type;
    typedef typename Multiplicable<first_value_type, second_value_type>::result_type result_value_type;
    
#if 0 // Just an idea
    typedef typename boost::mpl::if_<
                boost::mpl::and_<
	            boost::is_base_of<tag::sparse, typename traits::category<E1>::type>
	          , boost::is_base_of<tag::sparse, typename traits::category<E2>::type>
		> 
	      , compressed2D<result_value_type>
	      , dense2D<result_value_type, parameters<> >
            >::type                                       evaluated_result_type;
    
    // Convert into matrix

    operator evaluated_result_type() const
    {
	return evaluated_result_type(first * second);
    }
#endif

    mat_mat_times_expr( E1 const& v1, E2 const& v2 )
      : op_base( v1, v2 ), crtp_base(*this), first(v1), second(v2)
    {}

    // To prevent that cout << A * B prints the element-wise product, suggestion by Hui Li
    // It is rather inefficient, esp. for multiple products (complexity increases with the number of arguments :-!)
    //    or sparse matrices. 
    // Better compute your product first and print it then when compute time is an issue,
    // this is ONLY for convenience.
    result_value_type
    operator()(std::size_t r, std::size_t c) const
    {
	using math::zero;
	MTL_THROW_IF(num_cols(first) != num_rows(second), incompatible_size());

	result_value_type ref, sum(zero(ref));
	for (std::size_t i= 0; i < num_cols(first); i++)
	    sum+= first(r, i) * second(i, c);
	return sum;
    }

    
    result_value_type
    operator()(std::size_t r, std::size_t c)
    {
	return (*const_cast<const self*>(this))(r, c);
    }

    first_argument_type const&  first ;
    second_argument_type const& second ;
};

template <typename E1, typename E2>
std::size_t inline num_rows(const mat_mat_times_expr<E1, E2>& expr) 
{ return num_rows(expr.first); }

template <typename E1, typename E2>
std::size_t inline num_cols(const mat_mat_times_expr<E1, E2>& expr) 
{ return num_cols(expr.second); }

template <typename E1, typename E2>
std::size_t inline size(const mat_mat_times_expr<E1, E2>& expr) 
{ return num_rows(expr) * num_cols(expr); }

}} // Namespace mtl::matrix

#endif // MTL_MAT_MAT_TIMES_EXPR_INCLUDE
