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

#ifndef ITL_PC_IDENTITY_INCLUDE
#define ITL_PC_IDENTITY_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/itl/pc/solver.hpp>

namespace itl { namespace pc {

/// Identity preconditioner, i.e. no preconditioning and vector is just copied
/** Second template is just for a uniform interface with other preconditioners. **/
template <typename Matrix, typename Value= double>
class identity
{
  public:
    typedef typename mtl::Collection<Matrix>::value_type  value_type;
    typedef typename mtl::Collection<Matrix>::size_type   size_type;
    typedef identity                                      self;

    identity(const Matrix&) {}

    template <typename Vector>
    Vector solve(const Vector& x) const
    {
	mtl::vampir_trace<5032> tracer;
	return x;
    }

    template <typename VectorIn, typename VectorOut>
    void solve(const VectorIn& b, VectorOut& x) const
    {
	mtl::vampir_trace<5032> tracer;
	x= b;
    }

    template <typename Vector>
    Vector adjoint_solve(const Vector& x) const
    {
	mtl::vampir_trace<5034> tracer;
	return x;
    }

    template <typename VectorIn, typename VectorOut>
    void adjoint_solve(const VectorIn& b, VectorOut& x) const
    {
	mtl::vampir_trace<5034> tracer;
	x= b;
    }
}; 

template <typename Matrix, typename Vector>
solver<identity<Matrix>, Vector, false>
inline solve(const identity<Matrix>& P, const Vector& x)
{   return solver<identity<Matrix>, Vector, false>(P, x); }

template <typename Matrix, typename Vector>
solver<identity<Matrix>, Vector, true>
inline adjoint_solve(const identity<Matrix>& P, const Vector& x)
{   return solver<identity<Matrix>, Vector, true>(P, x); }


}} // namespace itl::pc

#endif // ITL_PC_IDENTITY_INCLUDE
