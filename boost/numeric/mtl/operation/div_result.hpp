/*
 *  div_result.h
 *  MTL
 *
 *  Created by Hui Li (huil@Princeton.EDU)
 *
 */

#ifndef MTL_DIV_RESULT_INCLUDE
#define MTL_DIV_RESULT_INCLUDE

#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>

namespace mtl { namespace traits {

template < typename Op1, typename Op2, typename DivOp > struct div_result_aux {};

/// Result type for dividing Op1 by Op2
/** Can be used in enable-if-style as type is only defined when appropriate **/
template < typename Op1, typename Op2 >
struct div_result 
  : public div_result_aux < Op1, Op2, typename ashape::div_op<typename ashape::ashape<Op1>::type, 
							      typename ashape::ashape<Op2>::type >::type > 
{};

/// Divide column vector by scalar
template < typename Op1, typename Op2 >
struct div_result_aux < Op1, Op2, ::mtl::ashape::cvec_scal_div >
{
    typedef typename vec::divide_by_view<Op1,Op2> type;
};

/// Divide row vector by scalar
template < typename Op1, typename Op2 >
struct div_result_aux < Op1, Op2, ::mtl::ashape::rvec_scal_div >
{
    typedef typename vec::divide_by_view<Op1,Op2> type;
};

/// Divide matrix by scalar
template < typename Op1, typename Op2 >
struct div_result_aux < Op1, Op2, ::mtl::ashape::mat_scal_div >
{
    typedef typename mat::divide_by_view<Op1,Op2> type;
};

}} // namespace mtl::traits


#endif
