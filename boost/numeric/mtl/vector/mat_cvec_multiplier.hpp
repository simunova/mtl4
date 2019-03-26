// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_VECTOR_MAT_CVEC_MULTIPLIER_INCLUDE
#define MTL_VECTOR_MAT_CVEC_MULTIPLIER_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/concept/std_concept.hpp>

#include <boost/numeric/mtl/vector/vec_expr.hpp>
#include <boost/numeric/mtl/vector/assigner.hpp>
#include <boost/numeric/mtl/vector/decrementer.hpp>
#include <boost/numeric/mtl/vector/incrementer.hpp>
#include <boost/numeric/mtl/operation/assign_mode.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace vec {

/// Helper class for delaying matrix vector multiplication, e.g. for matrix-free operators
template <typename Matrix, typename VectorIn>
struct mat_cvec_multiplier
  : assigner<mat_cvec_multiplier<Matrix, VectorIn> >,
    incrementer<mat_cvec_multiplier<Matrix, VectorIn> >,
    decrementer<mat_cvec_multiplier<Matrix, VectorIn> >,
    vec_expr<mat_cvec_multiplier<Matrix, VectorIn> >
{
    typedef typename Multiplicable<typename Collection<Matrix>::value_type,
				   typename Collection<VectorIn>::value_type>::result_type value_type;


    /// Construct with matrix \p A and vector \p v
    mat_cvec_multiplier(const Matrix& A, const VectorIn& v) : A(A), v(v) {}

    /// Assing the product to vector \p w, if possible directly within w's memory without copying
    template <typename VectorOut>
    void assign_to(VectorOut& w) const
    {
	vampir_trace<3068> tracer;
	A.mult(v, w, mtl::assign::assign_sum());
    }

    /// Increment vector \p w with the product, if possible directly within w's memory without copying
    template <typename VectorOut>
    void increment_it(VectorOut& w) const
    {
	vampir_trace<3068> tracer;
	A.mult(v, w, mtl::assign::plus_sum());
    }

    /// Decrement vector \p w with the product, if possible directly within w's memory without copying
    template <typename VectorOut>
    void decrement_it(VectorOut& w) const
    {
	vampir_trace<3068> tracer;
	A.mult(v, w, mtl::assign::minus_sum());
    }

#if 1
    /// Multiply vector \p w with the product, if possible directly within w's memory without copying
    template <typename VectorOut>
    void multiply_it(VectorOut& w) const
    {
    	vampir_trace<3068> tracer;
    	A.mult(v, w, mtl::assign::times_sum());
    }

    /// Divide vector \p w with the product, if possible directly within w's memory without copying
    template <typename VectorOut>
    void divide_it(VectorOut& w) const
    {
    	vampir_trace<3068> tracer;
    	A.mult(v, w, mtl::assign::divide_sum());
    }
#endif
    void delay_assign() const { }

    const Matrix&   A;
    const VectorIn& v;
};

template <typename Matrix, typename VectorIn>
inline std::size_t size(const mat_cvec_multiplier<Matrix, VectorIn>& m) 
{ return num_rows(m.A); }

template <typename Matrix, typename VectorIn>
inline std::size_t num_rows(const mat_cvec_multiplier<Matrix, VectorIn>& m) { return num_rows(m.A); }

template <typename Matrix, typename VectorIn>
inline std::size_t num_cols(const mat_cvec_multiplier<Matrix, VectorIn>&) { return 1; }
 
}} // namespace mtl::vector

#endif // MTL_VECTOR_MAT_CVEC_MULTIPLIER_INCLUDE
