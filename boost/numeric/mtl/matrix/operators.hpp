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

#ifndef MTL_MATRIX_OPERATORS_INCLUDE
#define MTL_MATRIX_OPERATORS_INCLUDE

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/and.hpp>

#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/is_multi_vector_expr.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/matrix/all_mat_expr.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>


namespace mtl { namespace mat {

template <typename E1, typename E2>
inline mat_mat_plus_expr<E1, E2>
operator+ (const mat_expr<E1>& e1, const mat_expr<E2>& e2)
{
    // do not add matrices with inconsistent value types
    MTL_STATIC_ASSERT((boost::is_same<typename ashape::ashape<E1>::type, 
				      typename ashape::ashape<E2>::type>::value), "Matrices have not consistent algebraic shape (i.e. nested types).");
    return mat_mat_plus_expr<E1, E2>(static_cast<const E1&>(e1), static_cast<const E2&>(e2));
}

// if enabled it has priority over previous functions because that performs upcast
template <typename E1, typename E2>
typename boost::enable_if_c<mtl::traits::is_multi_vector_expr<E1>::value && mtl::traits::is_multi_vector_expr<E2>::value, mv_mv_plus_expr<E1, E2> >::type
inline operator+(const E1& e1, const E2& e2)
{
    return mv_mv_plus_expr<E1, E2>(e1, e2);
}

// Specialization for multi_vector
// template <typename V1, typename V2>
// inline mv_mv_plus_expr<multi_vector<V1>, multi_vector<V2> >
// operator+(const multi_vector<V1>& m1, const multi_vector<V2>& m2)
// {
//     return mv_mv_plus_expr<multi_vector<V1>, multi_vector<V2> >(m1, m2);
// }

#if 0
// Planned for future optimizations on sums of dense matrix expressions
template <typename E1, typename E2>
inline dmat_dmat_plus_expr<E1, E2>
operator+ (const dmat_expr<E1>& e1, const dmat_expr<E2>& e2)
{
    // do not add matrices with inconsistent value types
    MTL_STATIC_ASSERT((boost::is_same<typename ashape::ashape<E1>::type, 
				      typename ashape::ashape<E2>::type>::value), "Matrices have not consistent algebraic shape (i.e. nested types).");
    return dmat_dmat_plus_expr<E1, E2>(static_cast<const E1&>(e1), static_cast<const E2&>(e2));
}
#endif


template <typename E1, typename E2>
inline mat_mat_minus_expr<E1, E2>
operator- (const mat_expr<E1>& e1, const mat_expr<E2>& e2)
{
    // do not add matrices with inconsistent value types
    MTL_STATIC_ASSERT((boost::is_same<typename ashape::ashape<E1>::type, 
				      typename ashape::ashape<E2>::type>::value), "Matrices have not consistent algebraic shape (i.e. nested types).");
    return mat_mat_minus_expr<E1, E2>(static_cast<const E1&>(e1), static_cast<const E2&>(e2));
}

// if enabled it has priority over previous functions because that performs upcast
template <typename E1, typename E2>
typename boost::enable_if_c<mtl::traits::is_multi_vector_expr<E1>::value && mtl::traits::is_multi_vector_expr<E2>::value, mv_mv_minus_expr<E1, E2> >::type
inline operator-(const E1& e1, const E2& e2)
{
    return mv_mv_minus_expr<E1, E2>(e1, e2);
}

template <typename E1, typename E2>
inline mat_mat_ele_times_expr<E1, E2>
ele_prod(const mat_expr<E1>& e1, const mat_expr<E2>& e2)
{
    // do not multiply matrices element-wise with inconsistent value types
    MTL_STATIC_ASSERT((boost::is_same<typename ashape::ashape<E1>::type, 
				      typename ashape::ashape<E2>::type>::value), "Matrices do not have consistent algebraic shape (i.e. nested types).");
    return mat_mat_ele_times_expr<E1, E2>(static_cast<const E1&>(e1), static_cast<const E2&>(e2));
}



}} // namespace mtl::matrix

#endif // MTL_MATRIX_OPERATORS_INCLUDE
