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

#ifndef ITL_BICGSTAB_2_INCLUDE
#define ITL_BICGSTAB_2_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/operation/resource.hpp>
#include <boost/numeric/mtl/operation/size.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/matrix/strict_upper.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/operation/orth.hpp>
#include <boost/numeric/mtl/operation/lazy.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#include <boost/numeric/itl/pc/identity.hpp>
#include <boost/numeric/itl/krylov/base_solver.hpp>

namespace itl {

/// Bi-Conjugate Gradient Stabilized(2)
template < typename LinearOperator, typename Vector, 
	   typename Preconditioner, typename Iteration >
int bicgstab_2(const LinearOperator &A, Vector &x, const Vector &b,
	       const Preconditioner &L, Iteration& iter)
{
    mtl::vampir_trace<7005> tracer;
    using mtl::size; using mtl::irange; using mtl::imax; using mtl::mat::strict_upper; using mtl::lazy;
    typedef typename mtl::Collection<Vector>::value_type Scalar;
    typedef typename mtl::Collection<Vector>::size_type  Size;

    if (size(b) == 0) throw mtl::logic_error("empty rhs vector");

    const size_t                l= 2;
    const Scalar                zero= math::zero(Scalar()), one= math::one(Scalar());
    Vector                      x0(resource(x)), y(resource(x));
    mtl::dense_vector<Vector>   r_hat(l+1,Vector(resource(x))), u_hat(l+1,Vector(resource(x)));

    // shift problem 
    x0= zero;
    r_hat[0]= b;
    if (two_norm(x) != zero) {
	r_hat[0]-= A * x;
	x0= x;
	x= zero;
    }

    Vector  r0_tilde(r_hat[0]/two_norm(r_hat[0]));
    y= solve(L, r_hat[0]);
    r_hat[0]= y;
    u_hat[0]= zero;

    Scalar                      rho_0(one), rho_1(zero), alpha(zero), Gamma(zero), beta(zero), omega(one); 
    mtl::mat::dense2D<Scalar>        tau(l+1, l+1);
    mtl::dense_vector<Scalar>   sigma(l+1), gamma(l+1), gamma_a(l+1), gamma_aa(l+1);

    while (! iter.finished(r_hat[0])) {
	++iter;
	rho_0= -omega * rho_0;

	for (Size j= 0; j < 2; ++j) {
	    rho_1= dot(r0_tilde, r_hat[j]); 
	    beta= alpha * rho_1/rho_0; rho_0= rho_1;

	    for (Size i= 0; i <= j; ++i)
		u_hat[i]= r_hat[i] - beta * u_hat[i];
      
	    y= A * u_hat[j];
	    u_hat[j+1]= solve(L, y);
	    Gamma= dot(r0_tilde, u_hat[j+1]); 
	    alpha= rho_0 / Gamma;

	    for (Size i= 0; i <= j; ++i)
		r_hat[i]-= alpha * u_hat[i+1];
      
	    if (iter.finished(r_hat[j])) {
		x+= x0;
		return iter;
	    }

	    y= A * r_hat[j]; 
	    r_hat[j+1]= solve(L, y);
	    x+= alpha * u_hat[0];
	}

	// mod GS (MR part)
	irange  i1m(1, imax);
	mtl::dense_vector<Vector>   r_hat_tail(r_hat[i1m]);
	tau[i1m][i1m]= orthogonalize_factors(r_hat_tail);
	for (Size j= 1; j <= l; ++j) 
	    gamma_a[j]= dot(r_hat[j], r_hat[0]) / tau[j][j];

	gamma[l]= gamma_a[l]; omega= gamma[l];
	if (omega == zero) return iter.fail(3, "bicg breakdown #2");

	// is this something like a tri-solve? 
	for (Size j= l-1; j > 0; --j) {
	    Scalar sum= zero;
	    for (Size i=j+1;i<=l;++i)
		sum += tau[j][i] * gamma[i];
	    gamma[j] = gamma_a[j] - sum;
	}

	gamma_aa[irange(1, l)]= strict_upper(tau[irange(1, l)][irange(1, l)]) * gamma[irange(2, l+1)] + gamma[irange(2, l+1)];

	x+= gamma[1] * r_hat[0];
	r_hat[0]-= gamma_a[l] * r_hat[l];
	u_hat[0]-= gamma[l] * u_hat[l];
	for (Size j=1; j < l; ++j) {
	    u_hat[0] -= gamma[j] * u_hat[j];
	    x+= gamma_aa[j] * r_hat[j];
	    r_hat[0] -= gamma_a[j] * r_hat[j];
	}
    }
    x+= x0; // convert to real solution and undo shift
    return iter;
}



/// Solver class for BiCGStab(2) method; right preconditioner ignored (prints warning if not identity)
/** Methods inherited from \ref base_solver. **/
template < typename LinearOperator, typename Preconditioner= pc::identity<LinearOperator>, 
	   typename RightPreconditioner= pc::identity<LinearOperator> >
class bicgstab_2_solver
  : public base_solver< bicgstab_2_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator >
{
    typedef base_solver< bicgstab_2_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator > base;
  public:
    /// Construct solver from a linear operator; generate (left) preconditioner from it
    explicit bicgstab_2_solver(const LinearOperator& A) : base(A), L(A) 
    {
	if (!pc::static_is_identity<RightPreconditioner>::value)
	    std::cerr << "Right Preconditioner ignored!" << std::endl;
    }

    /// Construct solver from a linear operator and (left) preconditioner
    bicgstab_2_solver(const LinearOperator& A, const Preconditioner& L) : base(A), L(L) 
    {
	if (!pc::static_is_identity<RightPreconditioner>::value)
	    std::cerr << "Right Preconditioner ignored!" << std::endl;
    }

    /// Solve linear system approximately as specified by \p iter
    template < typename HilbertSpaceX, typename HilbertSpaceB, typename Iteration >
    int solve(HilbertSpaceX& x, const HilbertSpaceB& b, Iteration& iter) const
    {
	return bicgstab_2(this->A, x, b, L, iter);
    }

  private:
    Preconditioner        L;
};



} // namespace itl

#endif // ITL_BICGSTAB_2_INCLUDE






