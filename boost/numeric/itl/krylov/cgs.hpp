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

#ifndef ITL_CGS_INCLUDE
#define ITL_CGS_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/resource.hpp>
#include <boost/numeric/mtl/operation/dot.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#include <boost/numeric/itl/krylov/base_solver.hpp>

namespace itl {

/// Conjugate Gradient Squared
template < typename LinearOperator, typename Vector, 
	   typename Preconditioner, typename Iteration >
int cgs(const LinearOperator &A, Vector &x, const Vector &b,
	const Preconditioner &M, Iteration& iter)
{
    mtl::vampir_trace<7007> tracer;
    typedef typename mtl::Collection<Vector>::value_type Scalar;
    Scalar     rho_1(0), rho_2(0), alpha(0), beta(0);
    Vector     p(resource(x)), phat(resource(x)), q(resource(x)), qhat(resource(x)), vhat(resource(x)),
	       u(resource(x)), uhat(resource(x)), r(b - A * x), rtilde= r;

    while (! iter.finished(r)) {
	++iter;
	rho_1= dot(rtilde, r);

	if (rho_1 == 0.) iter.fail(2, "cgs breakdown");

	if (iter.first())
	    p= u= r;
	else {
	    beta = rho_1 / rho_2;
	    u= r + beta * q;
	    p= u + beta * (q + beta * p);
	}

        vhat= A * Vector(solve(M, p));
	alpha = rho_1 / dot(rtilde, vhat);
	q= u - alpha * vhat;

	u+= q;
	uhat= solve(M, u);
	
	x+= alpha * uhat;
	qhat= A * uhat;
	r-= alpha * qhat;

	rho_2= rho_1;
    }
    return iter;
}

/// Solver class for CGS method; right preconditioner ignored (prints warning if not identity)
/** Methods inherited from \ref base_solver. **/
template < typename LinearOperator, typename Preconditioner= pc::identity<LinearOperator>, 
	   typename RightPreconditioner= pc::identity<LinearOperator> >
class cgs_solver
  : public base_solver< cgs_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator >
{
    typedef base_solver< cgs_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator > base;
  public:
    /// Construct solver from a linear operator; generate (left) preconditioner from it
    explicit cgs_solver(const LinearOperator& A) : base(A), L(A) 
    {
	if (!pc::static_is_identity<RightPreconditioner>::value)
	    std::cerr << "Right Preconditioner ignored!" << std::endl;
    }

    /// Construct solver from a linear operator and (left) preconditioner
    cgs_solver(const LinearOperator& A, const Preconditioner& L) : base(A), L(L) 
    {
	if (!pc::static_is_identity<RightPreconditioner>::value)
	    std::cerr << "Right Preconditioner ignored!" << std::endl;
    }

    /// Solve linear system approximately as specified by \p iter
    template < typename HilbertSpaceX, typename HilbertSpaceB, typename Iteration >
    int solve(HilbertSpaceX& x, const HilbertSpaceB& b, Iteration& iter) const
    {
	return cgs(this->A, x, b, L, iter);
    }

  private:
    Preconditioner        L;
};

} // namespace itl

#endif // ITL_CGS_INCLUDE



