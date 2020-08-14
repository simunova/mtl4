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

#ifndef ITL_BICG_INCLUDE
#define ITL_BICG_INCLUDE

#include <complex>
#include <boost/numeric/itl/itl_fwd.hpp>
#include <boost/numeric/itl/krylov/base_solver.hpp>

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/conj.hpp>
#include <boost/numeric/mtl/operation/resource.hpp>
#include <boost/numeric/mtl/operation/adjoint.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace itl {

/// Bi-Conjugate Gradient
template < typename LinearOperator, typename Vector, 
	   typename Preconditioner, typename Iteration >
int bicg(const LinearOperator &A, Vector &x, const Vector &b,
	 const Preconditioner &M, Iteration& iter)
{
    mtl::vampir_trace<7003> tracer;
    using mtl::conj;
    typedef typename mtl::Collection<Vector>::value_type Scalar;
    Scalar     rho_1(0), rho_2(0), alpha(0), beta(0);
    Vector     r(b - A * x), z(resource(x)), p(resource(x)), q(resource(x)),
 	       r_tilde(r), z_tilde(resource(x)), p_tilde(resource(x)), q_tilde(resource(x));

    while ( ! iter.finished(r)) {
	++iter;
	z= solve(M, r);
	z_tilde= adjoint_solve(M, r_tilde);
	rho_1= dot(z_tilde, z);

	if (rho_1 == 0.) return iter.fail(2, "bicg breakdown");
	if (iter.first()) {
	    p= z;
	    p_tilde= z_tilde;
	} else {
	    beta= rho_1 / rho_2;      
	    p= z + beta * p;
	    p_tilde= z_tilde + conj(beta) * p_tilde;
	}

	q= A * p;
	q_tilde= adjoint(A) * p_tilde;
	alpha= rho_1 / dot(p_tilde, q);

	x+= alpha * p;
	r-= alpha * q;
	r_tilde-= conj(alpha) * q_tilde;

	rho_2= rho_1;
    }
    return iter;
}

/// Solver class for BiCG method; right preconditioner ignored (prints warning if not identity)
/** Methods inherited from \ref base_solver. **/
template < typename LinearOperator, typename Preconditioner= pc::identity<LinearOperator>, 
	   typename RightPreconditioner= pc::identity<LinearOperator> >
class bicg_solver
  : public base_solver< bicg_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator >
{
    typedef base_solver< bicg_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator > base;
  public:
    /// Construct solver from a linear operator; generate (left) preconditioner from it
    explicit bicg_solver(const LinearOperator& A) : base(A), L(A) 
    {
	if (!pc::static_is_identity<RightPreconditioner>::value)
	    std::cerr << "Right Preconditioner ignored!" << std::endl;
    }

    /// Construct solver from a linear operator and (left) preconditioner
    bicg_solver(const LinearOperator& A, const Preconditioner& L) : base(A), L(L) 
    {
	if (!pc::static_is_identity<RightPreconditioner>::value)
	    std::cerr << "Right Preconditioner ignored!" << std::endl;
    }

    /// Solve linear system approximately as specified by \p iter
    template < typename HilbertSpaceX, typename HilbertSpaceB, typename Iteration >
    int solve(HilbertSpaceX& x, const HilbertSpaceB& b, Iteration& iter) const
    {
	return bicg(this->A, x, b, L, iter);
    }

  private:
    Preconditioner        L;
};

} // namespace itl

#endif // ITL_BICG_INCLUDE






