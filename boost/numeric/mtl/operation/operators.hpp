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

#ifndef MTL_OPERATORS_INCLUDE
#define MTL_OPERATORS_INCLUDE

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/and.hpp>

#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/matrix/operators.hpp>
//#include <boost/numeric/mtl/vector/operators.hpp>
#include <boost/numeric/mtl/operation/mult_result.hpp>
#include <boost/numeric/mtl/operation/div_result.hpp>
#include <boost/numeric/mtl/operation/dot.hpp>
#include <boost/numeric/mtl/operation/mat_cvec_times_expr.hpp>
#include <boost/numeric/mtl/matrix/all_mat_expr.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/true_copy.hpp>
#include <boost/numeric/mtl/vector/rvec_mat_times_expr.hpp>


namespace mtl { 

namespace mat {

    /// Multiplication for all supported types of operations
    /** Enable-if-like technique make sure that only called when properly defined **/
    template <typename Op1, typename Op2>
    typename mtl::traits::mult_result<typename mtl::traits::true_copy<Op1>::type, Op2>::type
    inline operator*(const Op1& op1, const Op2& op2)
    {
	
        return typename mtl::traits::mult_result<typename mtl::traits::true_copy<Op1>::type, Op2>::type(op1, op2);
    }



    /// Division of matrices and vectors by salars
    /** Enable-if-like technique make sure that only called when properly defined **/
    // added by Hui Li
    // enable_if_matrix shouldn't be needed but was nessecary in cuppen.hpp
    template < typename Op1, typename Op2 >
    typename mtl::traits::enable_if_matrix<Op1, typename mtl::traits::div_result<Op1,Op2>::type>::type
    inline operator/(const Op1& op1, const Op2& op2)
    {
        return typename mtl::traits::div_result<Op1,Op2>::type(op1,op2);
    }
	
} // namespace matrix


namespace vec {

    /// Multiplication for all supported types of operations
    /** Enable-if-like technique make sure that only called when properly defined **/
    template <typename Op1, typename Op2>
    typename mtl::traits::vec_mult_result<typename mtl::traits::true_copy<Op1>::type, Op2>::type
    inline operator*(const Op1& op1, const Op2& op2)
    {
        return typename mtl::traits::vec_mult_result<typename mtl::traits::true_copy<Op1>::type, Op2>::type(op1, op2);
    }

    /// Multiply row vector with column vector; result is scalar
    template <typename Op1, typename Op2>
    typename traits::lazy_enable_if_rvec_cvec_mult<Op1, Op2, detail::dot_result<Op1, Op2> >::type
    inline operator*(const Op1& op1, const Op2& op2)
    {
	return dot_real(op1, op2);
    }

    /// Division of matrices and vectors by salars
    /** Enable-if-like technique make sure that only called when properly defined **/
    // added by Hui Li
    template < typename Op1, typename Op2 >
    typename traits::div_result<Op1,Op2>::type
    inline operator/(const Op1& op1, const Op2& op2)
    {
        return typename traits::div_result<Op1,Op2>::type(op1,op2);
    }

    /// Compare two vectors for equality
    /** Enable-if makes sure that only called when properly defined **/
    template < typename Op1, typename Op2 >
    typename boost::enable_if<boost::mpl::and_<mtl::traits::is_vector<Op1>,
					       mtl::traits::is_vector<Op2> >, 
			      bool>::type 
    inline operator==(const Op1& op1, const Op2& op2)
    {
	std::size_t s1= num_rows(op1) * num_cols(op1), s2= num_rows(op2) * num_cols(op2); // ugly hack to fight with ADL
	if (s1 != s2)
	    return false;
	for (unsigned i= 0; i < s1; i++)
	    if (op1[i] != op2[i])
		return false;
	return true;
    }

    /// Compare two vectors for unequality
    /** Enable-if makes sure that only called when properly defined **/
    template < typename Op1, typename Op2 >
    typename boost::enable_if<boost::mpl::and_<mtl::traits::is_vector<Op1>,
					       mtl::traits::is_vector<Op2> >, 
			      bool>::type 
    inline operator!=(const Op1& op1, const Op2& op2)
    {
	return !(op1 == op2);
    }
	
} // namespace vector


} // namespace mtl

#endif // MTL_OPERATORS_INCLUDE
