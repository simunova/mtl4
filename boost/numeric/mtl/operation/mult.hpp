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

#ifndef MTL_MULT_INCLUDE
#define MTL_MULT_INCLUDE

#include <boost/numeric/mtl/config.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/flatcat.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/operation/dmat_dmat_mult.hpp>
#include <boost/numeric/mtl/operation/smat_smat_mult.hpp>
#include <boost/numeric/mtl/operation/smat_dmat_mult.hpp>
#include <boost/numeric/mtl/operation/mat_vec_mult.hpp>
#include <boost/numeric/mtl/operation/rvec_mat_mult.hpp> // Row vector times matrix
#include <boost/numeric/mtl/operation/mult_specialize.hpp>
#include <boost/numeric/mtl/operation/assign_mode.hpp>
#include <boost/numeric/mtl/operation/mult_assign_mode.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>

#include <boost/mpl/if.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace mat {


/// Multiplication: mult(a, b, c) computes c= a * b; 
/** The 3 types must be compatible, e.g. all three matrices or b and c are column vectors and a is a matrix.
    The dimensions are checked at compile time. **/
template <typename A, typename B, typename C>
typename mtl::traits::enable_if_matrix<A>::type
inline mult(const A& a, const B& b, C& c)
{
    vampir_trace<4010> tracer;
    MTL_DEBUG_THROW_IF(static_cast<const void*>(&a) == static_cast<const void*>(&c) 
		       || static_cast<const void*>(&b) == static_cast<const void*>(&c),
		       argument_result_conflict());

    // dispatch between matrices, vectors, and scalars
    using mtl::traits::shape_flatcat;
    gen_mult(a, b, c, assign::assign_sum(), shape_flatcat<A>(), shape_flatcat<B>(), shape_flatcat<C>());
                                            // typename category<A>::type(), typename category<B>::type(), typename category<C>::type());
}


/// Multiplication: mult_add(a, b, c) computes c+= a * b; 
/** The 3 types must be compatible, e.g. all three matrices or b and c are column vectors and a is a matrix.
    The dimensions are checked at compile time. **/
template <typename A, typename B, typename C>
typename mtl::traits::enable_if_matrix<A>::type
inline mult_add(const A& a, const B& b, C& c)
{
    vampir_trace<4010> tracer;
    // dispatch between matrices, vectors, and scalars
    using mtl::traits::shape_flatcat;
    gen_mult(a, b, c, assign::plus_sum(), shape_flatcat<A>(), shape_flatcat<B>(), shape_flatcat<C>());
                                          // typename category<A>::type(), typename category<B>::type(), typename category<C>::type());
}


/// Four term multiplication: mult(a, x, y, z) computes z= a * x + y; 
/** The 4 types must be compatible, i.e. a*x must be assignable to z and z must be incrementable by y.
    Right now, it is not more efficient than z= a * x; z+= y. For compatibility with MTL2. **/
template <typename A, typename X, typename Y, typename Z>
inline void mult(const A& a, const X& x, const Y& y, Z& z)
{
    vampir_trace<4010> tracer;
    mult(a, x, z);
    z+= y;
}


// Matrix multiplication
template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void gen_mult(const MatrixA& a, const MatrixB& b, MatrixC& c, Assign, tag::flat<tag::matrix>, tag::flat<tag::matrix>, tag::flat<tag::matrix>)
{
    vampir_trace<4011> tracer;
#if 1
    MTL_DEBUG_THROW_IF((const void*)&a == (const void*)&c || (const void*)&b == (const void*)&c, argument_result_conflict());
#else 
    if ((const void*)&a == (const void*)&c || (const void*)&b == (const void*)&c) {
	C tmp(num_rows(c), num_cols(c)); 
	mult(a, b, tmp);
	swap(C, tmp);
	return;
    }
#endif

    MTL_DEBUG_THROW_IF(num_rows(a) != num_rows(c) || num_cols(a) != num_rows(b) || num_cols(b) != num_cols(c), incompatible_size());
    // dispatch between dense and sparse
    using mtl::traits::sparsity_flatcat;
    mat_mat_mult(a, b, c, Assign(), sparsity_flatcat<MatrixA>(), sparsity_flatcat<MatrixB>(), sparsity_flatcat<MatrixC>());
		 // typename category<MatrixA>::type(), typename category<MatrixB>::type(), typename category<MatrixC>::type());
}


/// Dense matrix multiplication
/**  The function for dense matrix multiplication defines a default multiplication functor. 
     Alternatively the user can define his own functors for specific triplets of matrix types, 
     see detail::dmat_dmat_mult_specialize.
     The default functor for dense matrix multiplication is: 
     -# Use BLAS   if available, otherwise
     -# Recursive multiplication with:
        -# Platform optimized mult on blocks   if available, otherwise
        -# Tiled multiplication on blocks      if available, otherwise
        -# Naive multiplication on blocks
     -# Naive multiplication on entire matrices if recursion is not available
**/
template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void mat_mat_mult(const MatrixA& A, const MatrixB& b, MatrixC& c, Assign, tag::flat<tag::dense>, tag::flat<tag::dense>, tag::flat<tag::dense>)
{
    vampir_trace<4012> tracer;
    using assign::plus_sum; using assign::assign_sum; 

    static const unsigned long tiling1= detail::dmat_dmat_mult_tiling1<MatrixA, MatrixB, MatrixC>::value;
    static const unsigned long tiling2= detail::dmat_dmat_mult_tiling2<MatrixA, MatrixB, MatrixC>::value;
    typedef gen_tiling_dmat_dmat_mult_t<tiling1, tiling2, plus_sum>    tiling_mult_t;

    typedef gen_platform_dmat_dmat_mult_t<plus_sum, tiling_mult_t>     platform_mult_t;
    typedef gen_recursive_dmat_dmat_mult_t<platform_mult_t>            recursive_mult_t;
    typedef gen_blas_dmat_dmat_mult_t<assign_sum, recursive_mult_t>    blas_mult_t;
    typedef size_switch_dmat_dmat_mult_t<straight_dmat_dmat_mult_limit, tiling_mult_t, blas_mult_t>   variable_size_t;

    typedef fully_unroll_fixed_size_dmat_dmat_mult_t<Assign>           fully_unroll_t;
    typedef size_switch_dmat_dmat_mult_t<fully_unroll_dmat_dmat_mult_limit, fully_unroll_t, tiling_mult_t> fixed_size_t;

    static const bool all_static= mtl::traits::is_static<MatrixA>::value && mtl::traits::is_static<MatrixB>::value 
	                          && mtl::traits::is_static<MatrixC>::value;
    typedef static_switch_dmat_dmat_mult_t<all_static, fixed_size_t, variable_size_t>  default_functor_t;

    /// Use user-defined functor if provided (assign mode can be arbitrary)
    typedef typename boost::mpl::if_<
	detail::dmat_dmat_mult_specialize<MatrixA, MatrixB, MatrixC>
      , typename detail::dmat_dmat_mult_specialize<MatrixA, MatrixB, MatrixC>::type
      , default_functor_t
    >::type raw_functor_type;

    /// Finally substitute assign mode (consistently)
    typename assign::mult_assign_mode<raw_functor_type, Assign>::type functor;

    functor(A, b, c);
}

template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void mat_mat_mult(const MatrixA& A, const MatrixB& b, MatrixC& c, Assign, tag::flat<tag::dense>, tag::flat<tag::dense>, tag::flat<tag::sparse>)
{
    vampir_trace<4012> tracer;
    // This is a useless and extremely inefficient operation!!!!
    // We compute this with a dense matrix and copy the result back
    dense2D<typename Collection<MatrixC>::value_type, mat::parameters<> > c_copy(num_rows(c), num_cols(c));
    c_copy= c;
    mat_mat_mult(A, b, c_copy, Assign(), tag::flat<tag::dense>(), tag::flat<tag::dense>(), tag::flat<tag::dense>());
    c= c_copy;
}

/// Sparse matrix multiplication
template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void mat_mat_mult(const MatrixA& A, const MatrixB& b, MatrixC& c, Assign, tag::flat<tag::sparse>, tag::flat<tag::sparse>, tag::flat<tag::sparse>)
{
    vampir_trace<4012> tracer;
    smat_smat_mult(A, b, c, Assign(), typename OrientedCollection<MatrixA>::orientation(),
		   typename OrientedCollection<MatrixB>::orientation());
}

template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void mat_mat_mult(const MatrixA& A, const MatrixB& b, MatrixC& c, Assign, tag::flat<tag::sparse>, tag::flat<tag::sparse>, tag::flat<tag::dense>)
{
    vampir_trace<4012> tracer;
    // This is a useless and extremely inefficient operation!!!!
    // We compute this with a sparse matrix and copy the result back
    compressed2D<typename Collection<MatrixC>::value_type, mat::parameters<> > c_copy(num_rows(c), num_cols(c));
    c_copy= c;
    smat_smat_mult(A, b, c_copy, Assign(), typename OrientedCollection<MatrixA>::orientation(),
		   typename OrientedCollection<MatrixB>::orientation());
    c= c_copy;
}

/// Product of sparse times dense matrix
/**  This function (specialization of mult) is intended to multiply sparse matrices with multiple matrices
     gathered into a dense matrix.  Likewise, the resulting dense matrix corresponds to multiple vectors.
     The default functor for this operation is: 
     -# Use tiled multiplication      if available, otherwise
     -# Naive multiplication 
**/
template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void mat_mat_mult(const MatrixA& A, const MatrixB& b, MatrixC& c, Assign, tag::flat<tag::sparse>, tag::flat<tag::dense>, tag::flat<tag::dense>)
{
    vampir_trace<4012> tracer;
    using assign::plus_sum; using assign::assign_sum; 
    using namespace functor;

    // static const unsigned long tiling1= detail::dmat_dmat_mult_tiling1<MatrixA, MatrixB, MatrixC>::value;

    //typedef gen_smat_dmat_mult<Assign>                         default_functor_t;
    typedef gen_tiling_smat_dmat_mult<8, Assign>                         default_functor_t;

    // Finally substitute assign mode (consistently)
    // typename assign::mult_assign_mode<raw_functor_type, Assign>::type functor;

    default_functor_t functor;
    functor(A, b, c);
}

template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void mat_mat_mult(const MatrixA& A, const MatrixB& b, MatrixC& c, Assign, tag::flat<tag::sparse>, tag::flat<tag::dense>, tag::flat<tag::sparse>)
{
    vampir_trace<4012> tracer;
    // This is a useless and extremely inefficient operation!!!!
    // We compute this with a sparse matrix and copy the result back
    dense2D<typename Collection<MatrixC>::value_type, mat::parameters<> > c_copy(num_rows(c), num_cols(c));
    c_copy= c;
    mat_mat_mult(A, b, c_copy, Assign(), tag::flat<tag::sparse>(), tag::flat<tag::dense>(), tag::flat<tag::dense>());
    c= c_copy;
}


template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void mat_mat_mult(const MatrixA& A, const MatrixB& b, MatrixC& c, Assign, tag::flat<tag::dense>, tag::flat<tag::sparse>, tag::flat<tag::dense>)
{
    vampir_trace<4012> tracer;
    // This is could be a usefull operation, i.e. multiplying multiple row vectors with a sparse matrix
    // Might be supported in future
    // Now we compute this with a sparse matrix as first argument
    compressed2D<typename Collection<MatrixA>::value_type, mat::parameters<> > A_copy(num_rows(A), num_cols(A));
    A_copy= A;
    compressed2D<typename Collection<MatrixC>::value_type, mat::parameters<> > c_copy(num_rows(c), num_cols(c));
    c_copy= c;
    mat_mat_mult(A_copy, b, c_copy, Assign(), tag::flat<tag::sparse>(), tag::flat<tag::sparse>(), tag::flat<tag::sparse>());
    c= c_copy;
}



template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void mat_mat_mult(const MatrixA& A, const MatrixB& b, MatrixC& c, Assign, tag::flat<tag::dense>, tag::flat<tag::sparse>, tag::flat<tag::sparse>)
{
    vampir_trace<4012> tracer;
    // This is not a usefull operation, because the result is dense
    // Now we compute this with a sparse matrix as first argument
    compressed2D<typename Collection<MatrixA>::value_type, mat::parameters<> > A_copy(num_rows(A), num_cols(A));
    A_copy= A;
    mat_mat_mult(A_copy, b, c, Assign(), tag::flat<tag::sparse>(), tag::flat<tag::sparse>(), tag::flat<tag::sparse>());
}



// Matrix vector multiplication
template <typename Matrix, typename VectorIn, typename VectorOut, typename Assign>
inline void gen_mult(const Matrix& A, const VectorIn& v, VectorOut& w, Assign, tag::flat<tag::matrix>, tag::flat<tag::col_vector>, tag::flat<tag::col_vector>)
{
    vampir_trace<4011> tracer;
    // Vector must be column vector
    // If vector is row vector then matrix must have one column and the operation is a outer product
    //   -> result should be a matrix too

    // Check if element types are compatible (in contrast to tag dispatching, nesting is considered here)
    MTL_STATIC_ASSERT((boost::is_same< typename ashape::mult_op<typename ashape::ashape<Matrix>::type, 
			                                          typename ashape::ashape<VectorIn>::type >::type,
			                 ::mtl::ashape::mat_cvec_mult
				     >::value),
		      "The type nesting of the arguments does not allow for a consistent matrix vector product.");


#if 1
    MTL_DEBUG_THROW_IF((void*)&v == (void*)&w, argument_result_conflict());
#else
    if ((void*)&v == (void*)&w) {
	VectorOut tmp(size(w)); 
	mult(A, b, tmp);
	swap(w, tmp);
	return;
    }
#endif
    // w.checked_change_dim(num_rows(A)); // destroys distribution in parallel -> dimension changed in assignment
    MTL_DEBUG_THROW_IF(num_rows(A) != mtl::size(w), incompatible_size());
    MTL_DEBUG_THROW_IF(num_cols(A) != mtl::size(v), incompatible_size());

    mat_cvec_mult(A, v, w, Assign(), mtl::traits::mat_cvec_flatcat<Matrix>());
}


// Vector matrix multiplication
template <typename VectorIn, typename Matrix, typename VectorOut, typename Assign>
inline void gen_mult(const VectorIn& v, const Matrix& A, VectorOut& w, Assign, tag::flat<tag::row_vector>, tag::flat<tag::matrix>, tag::flat<tag::row_vector>)
{
    vampir_trace<4011> tracer;
    // Vector must be column vector
    // If vector is row vector then matrix must have one column and the operation is a outer product
    //   -> result should be a matrix too

    // Check if element types are compatible (in contrast to tag dispatching, nesting is considered here)
    MTL_STATIC_ASSERT((boost::is_same< typename ashape::mult_op<typename ashape::ashape<VectorIn>::type, 
			                                          typename ashape::ashape<Matrix>::type >::type,
			                 ::mtl::ashape::rvec_mat_mult
			             >::value),
		      "The type nesting of the arguments does not allow for a consistent matrix vector product.");

#if 1
    MTL_DEBUG_THROW_IF((void*)&v == (void*)&w, argument_result_conflict());
#else
    if ((void*)&v == (void*)&w) {
	VectorOut tmp(size(w)); 
	mult(A, b, tmp);
	swap(w, tmp);
	return;
    }
#endif
    // w.checked_change_dim(num_cols(A));
    w.checked_change_resource(v);
    MTL_DEBUG_THROW_IF(num_cols(v) != num_rows(A), incompatible_size());

    // same dispatching criterion as mat_cvec_mult (until now)
    rvec_mat_mult(v, A, w, Assign(), typename mtl::traits::mat_cvec_flatcat<Matrix>::type());
}





}} // namespace mtl::matrix

#endif // MTL_MULT_INCLUDE
