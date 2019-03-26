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

#ifndef MTL_ASHAPE_INCLUDE
#define MTL_ASHAPE_INCLUDE

#include <vector>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/root.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>

// Not elegant but necessary to treat ITL types right
#include <boost/numeric/itl/itl_fwd.hpp>

#ifdef MTL_WITH_INITLIST
# include <initializer_list>
#endif

namespace mtl { 

/// Namespace for algebraic shapes; used for sophisticated dispatching between operations
namespace ashape {

// forward declaration
template <typename T> struct ashape_aux;

/// Tag for arbitrary algebraic shape
struct universe {};

// Types (tags)
/// Scalar algebraic shape
struct scal : universe {};

/// Non-scalar algebraic shape
struct nonscal : universe {};
/// Row vector as algebraic shape
template <typename Value> struct rvec : nonscal {};
/// Column vector as algebraic shape
template <typename Value> struct cvec : nonscal {};
/// Matrix as algebraic shape
template <typename Value> struct mat : nonscal {};
/// Undefined shape, e.g., for undefined results of operations
struct ndef {};
/// Future shape, i.e. after appropriate evaluation it will have the shape \p Value
template <typename Value> struct future : nonscal {};

/// Meta-function for algebraic shape of T
/** Unknown types are treated like scalars. ashape of collections are template
    parameterized with ashape of their elements, e.g., ashape< matrix < vector < double > > >::type is
    mat< rvec < scal > > >. 
    Implemented with ashape_aux after type is cleaned up with mtl::traits::root.
**/
template <typename T>
struct ashape
  : ashape_aux<typename mtl::traits::root<T>::type> {};

template <typename T>
struct ashape_aux
{
    typedef scal type;
};

/// Vectors must be distinguished between row and column vectors
template <typename Value, typename Parameters>
struct ashape_aux<mtl::vec::dense_vector<Value, Parameters> >
{
    typedef typename boost::mpl::if_<
  boost::is_same<typename Parameters::orientation, row_major>
      , rvec<typename ashape<Value>::type>
      , cvec<typename ashape<Value>::type>
    >::type type;
};

/// Same as dense vector
template <typename Value, typename Parameters>
struct ashape_aux<vec::strided_vector_ref<Value, Parameters> >
  : ashape<vec::dense_vector<Value, Parameters> > {};

/// Same as dense vector
template <typename Value, typename Parameters>
struct ashape_aux<vec::sparse_vector<Value, Parameters> >
  : ashape<vec::dense_vector<Value, Parameters> > {};

/// One-dimensional arrays have rvec ashape; 2D arrays are matrices see below
template <typename Value, unsigned Rows>
struct ashape_aux<Value[Rows]>
{
    typedef rvec<typename ashape<Value>::type> type;
};
   
#ifdef MTL_WITH_INITLIST
/// Non-nested initializer_list have rvec ashape, nested lists are matrices see below
template <typename Value>
struct ashape_aux<std::initializer_list<Value> >
{
    typedef rvec<typename ashape<Value>::type> type;
};
#endif

/// std::vectors have rvec ashape
template <typename Value, typename Allocator>
struct ashape_aux<std::vector<Value, Allocator> >
{
    typedef rvec<typename ashape<Value>::type> type;
};
   
/// One-dimensional arrays have rvec ashape; 2D arrays are matrices see below
template <typename Value>
struct ashape_aux<Value*>
{
    typedef rvec<typename ashape<Value>::type> type;
};
   
template <typename E1, typename E2, typename SFunctor>
struct ashape_aux< vec::vec_vec_pmop_expr<E1, E2, SFunctor> >
{
    MTL_STATIC_ASSERT((boost::is_same<typename ashape<E1>::type, 
                      typename ashape<E2>::type>::value), "Operands must have same algebraic shape.");
    typedef typename ashape<E1>::type type;
};

template <typename E1, typename E2, typename SFunctor>
struct ashape_aux< vec::vec_vec_op_expr<E1, E2, SFunctor> >
{
#if 0 // not sure if this is true in all operations
    MTL_STATIC_ASSERT((boost::is_same<typename ashape<E1>::type, 
              typename ashape<E2>::type>::value), "Operands must have same algebraic shape.");
#endif
    typedef typename ashape<E1>::type type;
};

template <typename E1, typename E2, typename SFunctor>
struct ashape_aux< vec::vec_vec_aop_expr<E1, E2, SFunctor> >
{
    typedef typename ashape<E1>::type type;
};

template <typename E1, typename E2, typename SFunctor>
struct ashape_aux< vec::vec_scal_aop_expr<E1, E2, SFunctor> >
{
    typedef typename ashape<E1>::type type;
};

template <typename Vector>
struct ashape_aux< vec::vec_const_ref_expr<Vector> >
{
    typedef typename ashape<Vector>::type type;
};


// ========
// Matrices
// ========

template <typename Value, typename Parameters>
struct ashape_aux<mtl::mat::compressed2D<Value, Parameters> >
{
    typedef mat<typename ashape<Value>::type> type;
};

template <typename Value, typename Parameters>
struct ashape_aux<mtl::mat::coordinate2D<Value, Parameters> >
{
    typedef mat<typename ashape<Value>::type> type;
};

template <typename Value, typename Parameters>
struct ashape_aux<mtl::mat::sparse_banded<Value, Parameters> >
{
    typedef mat<typename ashape<Value>::type> type;
};

template <typename Value, typename Parameters>
struct ashape_aux<mtl::mat::dense2D<Value, Parameters> >
{
    typedef mat<typename ashape<Value>::type> type;
};
   
template <typename Value, std::size_t Mask, typename Parameters>
struct ashape_aux<mtl::mat::morton_dense<Value, Mask, Parameters> >
{
    typedef mat<typename ashape<Value>::type> type;
};

template <typename Functor>
struct ashape_aux<mtl::mat::implicit_dense<Functor> >
{
    typedef mat<typename ashape<typename Functor::result_type>::type> type;
};

/// Two-dimensional arrays have mat ashape; 1D arrays are vectors see above
template <typename Value, unsigned Rows, unsigned Cols>
struct ashape_aux<Value[Rows][Cols]>
{
    typedef mat<typename ashape<Value>::type> type;
};

/// Two-dimensional arrays have mat ashape; 1D arrays are vectors see above
template <typename Value, unsigned Cols>
struct ashape_aux<Value (*)[Cols]>
{
    typedef mat<typename ashape<Value>::type> type;
};

#ifdef MTL_WITH_INITLIST
/// Nested initializer_list are matrices, non-nested are vectors see above
template <typename Value>
struct ashape_aux<std::initializer_list<std::initializer_list<Value> > >
{
    typedef mat<typename ashape<Value>::type> type;
};
#endif

template <typename Vector>
struct ashape_aux<mtl::mat::multi_vector<Vector> >
{
    typedef mat<typename ashape<typename mtl::Collection<mtl::mat::multi_vector<Vector> >::value_type>::type> type;
};
  
template <typename Value>
struct ashape_aux<mtl::mat::element_structure<Value> >
{
   typedef mat<typename ashape<Value>::type> type;
};

template <typename Value, typename Parameters>
struct ashape_aux<mtl::mat::ell_matrix<Value, Parameters> >
{
   typedef mat<typename ashape<Value>::type> type;
};

 
template <typename Vector>
struct ashape_aux<mtl::mat::multi_vector_range<Vector> >
{
    typedef mat<typename ashape<typename mtl::Collection<mtl::mat::multi_vector_range<Vector> >::value_type>::type> type;
};

template <> struct ashape_aux<mtl::mat::identity2D> 
{  
    typedef nonscal type;    
};

template <typename E1, typename E2, typename SFunctor>
struct ashape_aux< mtl::mat::mat_mat_op_expr<E1, E2, SFunctor> >
{
    MTL_STATIC_ASSERT((boost::is_same<typename ashape<E1>::type, 
              typename ashape<E2>::type>::value), "Operands must have same algebraic shape.");
    typedef typename ashape<E1>::type type;
};

template <typename Vector1, typename Vector2>
struct ashape< mtl::mat::outer_product_matrix<Vector1, Vector2> >
{
    // BOOST_STATIC_ASSERT((boost::is_same<typename ashape<E1>::type, 
    //                       typename transposed_shape<typename ashape<E2>::type>::type>::value));    
    typedef mat<typename ashape<typename mtl::Collection<Vector1>::value_type>::type> type;
};

template <typename Matrix, typename VectorIn> 
struct ashape< mtl::vec::mat_cvec_multiplier<Matrix, VectorIn> >
{
    typedef cvec<scal> type;
};

// =====
// Views
// =====

template <typename Functor, typename Coll>
struct ashape_aux<mtl::mat::map_view<Functor, Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Functor, typename Coll>
struct ashape_aux<mtl::vec::map_view<Functor, Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Vector, typename Exponent>
struct ashape_aux<mtl::vec::pow_by_view<Vector, Exponent> >
{
    typedef typename ashape<Vector>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::conj_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::real_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::imag_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::negate_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::abs_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::acos_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::acosh_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::asin_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::asinh_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::atan_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::atanh_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::cos_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::cosh_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::sin_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::sinh_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::tan_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::tanh_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::ceil_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::floor_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::log_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::log10_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::exp_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::exp10_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::sqrt_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::rsqrt_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::signum_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};


# ifdef MTL_WITH_MATH_ELEVEN    
template <typename Coll>
struct ashape_aux<mtl::vec::round_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::trunc_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::log2_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::exp2_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::erf_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};

template <typename Coll>
struct ashape_aux<mtl::vec::erfc_view<Coll> >
{
    typedef typename ashape<Coll>::type type;
};
# endif




#if 1
// shouldn't be needed 
template <typename Coll>
struct ashape_aux<mtl::mat::transposed_view<const mtl::mat::conj_view<Coll> > >
{
    typedef typename ashape<Coll>::type type;
};
#endif

template <typename Matrix>
struct ashape_aux<mtl::mat::transposed_view<Matrix> >
{
    typedef typename ashape<Matrix>::type type;
};

template <typename Matrix>
struct ashape_aux<mtl::mat::banded_view<Matrix> >
{
    typedef typename ashape<Matrix>::type type;
};

template <typename Matrix>
struct ashape_aux<mtl::mat::indirect<Matrix> >
{
    typedef typename ashape<Matrix>::type type;
};

// Rule out other types as algebraic shape
template <typename IFStream, typename OFStream>
struct ashape_aux<io::matrix_file<IFStream, OFStream> > 
{
    typedef ndef type;
};


// =====================
// Shapes of products:
// =====================

// a) The result's shape
// b) Classify operation in terms of shape

// Operation types:

struct scal_scal_mult {};
struct cvec_rvec_mult {}; // outer product
struct rvec_cvec_mult {}; // inner product (without conj)
struct rvec_rvec_mult {}; // element-wise product
struct cvec_cvec_mult {}; // element-wise product
struct rvec_mat_mult {};
struct mat_cvec_mult {};
struct mat_mat_mult {};
struct scal_rvec_mult {};
struct scal_cvec_mult {};
struct scal_mat_mult {};
struct rvec_scal_mult {};
struct cvec_scal_mult {};
struct mat_scal_mult {};



// =====================
// Results of operations
// =====================

/* 
      s  cv  rv   m
-------------------
 s |  s  cv* rv*  m*
cv | cv*  x   m   x
rv | rv*  s   x  rv
 m |  m* cv   x   m 

 * only on outer level, forbidden for elements of collections

*/

// Results for elements of collections, i.e. scalar * matrix (vector) are excluded


/// Algebraic shape of multiplication's result when elements of collections are multiplied.
/** The types are the same as for multiplications of entire collections except that scalar *
    matrix (or vector) is excluded to avoid ambiguities. 
    emult_shape <Shape1, Shape2> is only properly defined if emult_op <Shape1, Shape2>::type is not ndef!
**/
template <typename Shape1, typename Shape2>
struct emult_shape
{
    typedef ndef type;
};

/// Type of operation when values of Shape1 and Shape2 are multiplied (so far only for elements of collections)
/** The types are the same as for multiplications of entire collections except that scalar *
    matrix (or vector) is excluded to avoid ambiguities. **/
template <typename Shape1, typename Shape2>
struct emult_op
{
    typedef ndef type;
};


// Scalar * scalar -> scalar
template <>
struct emult_shape<scal, scal>
{
    typedef scal type;
};

template <>
struct emult_op<scal, scal>
{
    typedef scal_scal_mult type;
};

// Column times row vector, i.e. outer product
template <typename Value1, typename Value2>
struct emult_shape<cvec<Value1>, rvec<Value2> >
{
    typedef mat<typename emult_shape<Value1, Value2>::type> type;
};

template <typename Value1, typename Value2>
struct emult_op<cvec<Value1>, rvec<Value2> >
{
    // if product of elements is undefined then product is undefined too
    typedef typename boost::mpl::if_<
  boost::is_same<typename emult_op<Value1, Value2>::type, ndef>
      , ndef
      , cvec_rvec_mult
    >::type type;
};


// Row times column vector, i.e. inner product (without conj)
template <typename Value1, typename Value2>
struct emult_shape<rvec<Value1>, cvec<Value2> >
{
    typedef typename emult_shape<Value1, Value2>::type type;
};

template <typename Value1, typename Value2>
struct emult_op<rvec<Value1>, cvec<Value2> >
{
    // if product of elements is undefined then product is undefined too
    typedef typename boost::mpl::if_<
  boost::is_same<typename emult_op<Value1, Value2>::type, ndef>
      , ndef
      , rvec_cvec_mult
    >::type type;
};

#ifndef MTL_WITHOUT_VECTOR_ELE_OPS
template <typename Value1, typename Value2>
struct emult_op<rvec<Value1>, rvec<Value2> >
{
    // if product of elements is undefined then product is undefined too
    typedef typename boost::mpl::if_<
  boost::is_same<typename emult_op<Value1, Value2>::type, ndef>
      , ndef
      , rvec_rvec_mult
    >::type type;
};

template <typename Value1, typename Value2>
struct emult_op<cvec<Value1>, cvec<Value2> >
{
    // if product of elements is undefined then product is undefined too
    typedef typename boost::mpl::if_<
  boost::is_same<typename emult_op<Value1, Value2>::type, ndef>
      , ndef
      , cvec_cvec_mult
    >::type type;
};
#endif // MTL_WITHOUT_VECTOR_ELE_OPS


// Row vector times matrix
template <typename Value1, typename Value2>
struct emult_shape<rvec<Value1>, mat<Value2> >
{
    typedef rvec<typename emult_shape<Value1, Value2>::type> type;
};


template <typename Value1, typename Value2>
struct emult_op<rvec<Value1>, mat<Value2> >
{
    // if product of elements is undefined then product is undefined too
    typedef typename boost::mpl::if_<
  boost::is_same<typename emult_op<Value1, Value2>::type, ndef>
      , ndef
      , rvec_mat_mult
    >::type type;
};

// Matrix times column vector
template <typename Value1, typename Value2>
struct emult_shape<mat<Value1>, cvec<Value2> >
{
    typedef cvec<typename emult_shape<Value1, Value2>::type> type;
};

template <typename Value1, typename Value2>
struct emult_op<mat<Value1>, cvec<Value2> >
{
    // if product of elements is undefined then product is undefined too
    typedef typename boost::mpl::if_<
  boost::is_same<typename emult_op<Value1, Value2>::type, ndef>
      , ndef
      , mat_cvec_mult
    >::type type;
};


// Matrix product
template <typename Value1, typename Value2>
struct emult_shape<mat<Value1>, mat<Value2> >
{
    typedef mat<typename emult_shape<Value1, Value2>::type> type;
};

template <typename Value1, typename Value2>
struct emult_op<mat<Value1>, mat<Value2> >
{
    // if product of elements is undefined then product is undefined too
    typedef typename boost::mpl::if_<
  boost::is_same<typename emult_op<Value1, Value2>::type, ndef>
      , ndef
      , mat_mat_mult
    >::type type;
};


// Results for entire collections, i.e. scalar * matrix (vector) are allowed

// Multiplying collections as emult

template  <typename Shape1, typename Shape2>
struct mult_shape
    : public emult_shape<Shape1, Shape2>
{};

template  <typename Shape1, typename Shape2>
struct mult_op
    : public emult_op<Shape1, Shape2>
{};

// Scale collection from left

template <typename Shape2>
struct mult_shape<scal, Shape2>
{
    typedef Shape2 type;
};

template <typename Value2>
struct mult_op<scal, rvec<Value2> >
{
    typedef scal_rvec_mult type;
};

template <typename Value2>
struct mult_op<scal, cvec<Value2> >
{
    typedef scal_cvec_mult type;
};

template <typename Value2>
struct mult_op<scal, mat<Value2> >
{
    typedef scal_mat_mult type;
};

// Scale collection from right

template <typename Shape1>
struct mult_shape<Shape1, scal>
{
    typedef Shape1 type;
};

template <typename Value1>
struct mult_op<rvec<Value1>, scal>
{
    typedef rvec_scal_mult type;
};

template <typename Value1>
struct mult_op<cvec<Value1>, scal>
{
    typedef cvec_scal_mult type;
};

template <typename Value1>
struct mult_op<mat<Value1>, scal>
{
    typedef mat_scal_mult type;
};

// Arbitration
template <>
struct mult_shape<scal, scal>
{
    typedef scal type;
};


// Needs to be verified for nested matrix types, cf. #140
template <typename E1, typename E2>
struct ashape< mtl::mat::mat_mat_times_expr<E1, E2> >
{
    // typedef typename ashape<E1>::type type;
    typedef typename mult_shape<typename ashape<E1>::type, 
        typename ashape<E2>::type>::type type;
};


template <typename E1, typename E2>
struct ashape< mat_cvec_times_expr<E1, E2> >
{
    // Resulting vector has the same shape as the multiplied
    typedef typename ashape<E2>::type type;
};

template <typename E1, typename E2>
struct ashape< vec::rvec_mat_times_expr<E1, E2> >
{
    // Resulting vector has the same shape as the multiplied
    typedef typename ashape<E1>::type type;
};


// added by Hui Li (below) -----------------------------------------

// =====================
// Shapes of divisions:
// =====================

// Operation types:

struct scal_scal_div {};
struct cvec_scal_div {};
struct rvec_scal_div {};
struct mat_scal_div {};

template < typename Shape1, typename Shape2 >
struct div_shape
{
  typedef ndef type;
};

template < typename Shape1, typename Shape2 >
struct div_op
{
  typedef ndef type;
};

template <>
struct div_shape<scal,scal>
{
  typedef scal type;
};

template<>
struct div_op<scal,scal>
{
  typedef scal type;
};

template < typename Value1 >
struct div_shape < rvec<Value1>, scal >
{
  typedef typename boost::mpl::if_<
    typename boost::is_same<typename div_shape<Value1,scal>::type,ndef>::type,
    ndef,
    rvec<typename div_shape<Value1,scal>::type>
  >::type type;
};

template < typename Value1 >
struct div_op< rvec<Value1>, scal >
{
  typedef typename boost::mpl::if_<
    typename boost::is_same<typename div_shape<rvec<Value1>,scal>::type,ndef>::type,
    ndef,
    rvec_scal_div
  >::type type;
};

template < typename Value1 >
struct div_shape < cvec<Value1>, scal >
{
  typedef typename boost::mpl::if_<
    typename boost::is_same<typename div_shape<Value1,scal>::type,ndef>::type,
    ndef,
    cvec<typename div_shape<Value1,scal>::type>
  >::type type;
};

template < typename Value1 >
struct div_op< cvec<Value1>, scal >
{
  typedef typename boost::mpl::if_<
    typename boost::is_same<typename div_shape<cvec<Value1>,scal>::type,ndef>::type,
    ndef,
    cvec_scal_div
  >::type type;
};

template < typename Value1 >
struct div_shape < mat<Value1>, scal >
{
  typedef typename boost::mpl::if_<
    typename boost::is_same<typename div_shape<Value1,scal>::type,ndef>::type,
    ndef,
    mat<typename div_shape<Value1,scal>::type>
  >::type type;
};

template < typename Value1 >
struct div_op < mat<Value1>, scal >
{
  typedef typename boost::mpl::if_<
    typename boost::is_same<typename div_shape<mat<Value1>,scal>::type,ndef>::type,
    ndef,
    mat_scal_div
  >::type type;
};
  
// added by Hui Li (above) -----------------------------------------

// ==================== ITL types ==================================

template <typename PC, typename Vector, bool Adjoint>
struct ashape<itl::pc::solver<PC, Vector, Adjoint> >
{
    typedef future<cvec<scal> > type; // might be a problem with nested matrices and vectors
};


}} // namespace mtl::ashape

#endif // MTL_ASHAPE_INCLUDE
