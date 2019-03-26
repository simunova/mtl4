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

// Written by Cornelius Steinhardt


#ifndef ITL_TFQMR_INCLUDE
#define ITL_TFQMR_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/operation/resource.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#include <boost/numeric/itl/krylov/base_solver.hpp>

namespace itl {

/// Transposed-free Quasi-minimal residual
template < typename Matrix, typename Vector,
	   typename LeftPreconditioner, typename RightPreconditioner, typename Iteration >
int tfqmr(const Matrix &A, Vector &x, const Vector &b, const LeftPreconditioner &L, 
	  const RightPreconditioner &R, Iteration& iter)
{
    mtl::vampir_trace<7009> tracer;
    using math::reciprocal; using mtl::size;
    typedef typename mtl::Collection<Vector>::value_type Scalar;

    if (size(b) == 0) throw mtl::logic_error("empty rhs vector");

    const Scalar                zero= math::zero(Scalar()), one= math::one(Scalar());
    Scalar                      theta(zero), eta(zero), tau, rho, rhon, sigma, alpha, beta, c;
    Vector                      rt(b - A*Vector(solve(R, x))) /* shift x= R*x */, r(solve(L, rt)), u1(resource(x)), u2(resource(x)), 
                                y1(resource(x)), y2(resource(x)), w(resource(x)), d(resource(x), zero), v(resource(x));

    if (iter.finished(rt))
	return iter;
    y1= w= r;
    rt= A * Vector(solve(R, y1));
    u1= v= solve(L,rt);
    tau= two_norm(r);
    rho= tau * tau;

    // TFQMR iteration
    while (! iter.finished(tau)) {
	++iter;
	sigma= dot(r,v);
        if (sigma == zero)
	    return iter.fail(1, "tfgmr breakdown, sigma=0 #1");
        alpha= rho / sigma;

        // inner loop
        for(int j=1; j < 3; j++) {
            if (j == 1) {
                w-= alpha * u1;
                d= y1+ (theta * theta * eta / alpha) * d;
	    } else {
                y2= y1 - alpha * v;
                rt= A * Vector(solve(R, y2));
                u2= solve(L, rt);
                w-= alpha * u2;
                d= y2 + (theta * theta * eta / alpha) * d;
            }
            theta= two_norm(w) / tau;
            c= reciprocal(sqrt(one + theta*theta));
            tau*= theta * c;
            eta= c * c * alpha;
            x+= eta * d;
        } // end inner loop
        if (rho == zero)
            return iter.fail(1, "tfgmr breakdown, rho=0 #2");
        rhon= dot(r,w);
        beta= rhon/rho;
        rho= rhon;
        y1= w + beta*y2;
        rt= A * Vector(solve(R, y1));
        u1= solve(L, rt);
        v= u1 + beta*(u2 + beta*v);
        rt= A * x - b;
    }
    //shift back
    x= solve(R, x);
    return iter;
}

/// Solver class for Transposed-free quasi-minimal residual method; right preconditioner ignored (prints warning if not identity)
/** Methods inherited from \ref base_solver. **/
template < typename LinearOperator, typename Preconditioner= pc::identity<LinearOperator>, 
	   typename RightPreconditioner= pc::identity<LinearOperator> >
class tfqmr_solver
  : public base_solver< tfqmr_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator >
{
    typedef base_solver< tfqmr_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator > base;
  public:
    /// Construct solver from a linear operator; generate (left) preconditioner from it
    explicit tfqmr_solver(const LinearOperator& A) : base(A), L(A), R(A) {}

    /// Construct solver from a linear operator and left preconditioner
    tfqmr_solver(const LinearOperator& A, const Preconditioner& L) : base(A), L(L), R(A) {}

    /// Construct solver from a linear operator and left preconditioner
    tfqmr_solver(const LinearOperator& A, const Preconditioner& L, const RightPreconditioner& R) 
      : base(A), L(L), R(R) {}

    /// Solve linear system approximately as specified by \p iter
    template < typename HilbertSpaceX, typename HilbertSpaceB, typename Iteration >
    int solve(HilbertSpaceX& x, const HilbertSpaceB& b, Iteration& iter) const
    {
	return tfqmr(this->A, x, b, L, R, iter);
    }

  private:
    Preconditioner        L;
    RightPreconditioner   R;
};

} // namespace itl

#endif // ITL_TFQMR_INCLUDE
