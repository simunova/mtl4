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

#ifndef ITL_BICGSTAB_INCLUDE
#define ITL_BICGSTAB_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/operation/resource.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#include <boost/numeric/itl/utility/exception.hpp>
#include <boost/numeric/itl/krylov/base_solver.hpp>

namespace itl {

///  Bi-Conjugate Gradient Stabilized
template < class LinearOperator, class HilbertSpaceX, class HilbertSpaceB, 
	   class Preconditioner, class Iteration >
int bicgstab(const LinearOperator& A, HilbertSpaceX& x, const HilbertSpaceB& b, 
	     const Preconditioner& M, Iteration& iter)
{
  typedef typename mtl::Collection<HilbertSpaceX>::value_type Scalar;
  typedef HilbertSpaceX                                       Vector;
  mtl::vampir_trace<7004> tracer;

  Scalar     rho_1(0), rho_2(0), alpha(0), beta(0), gamma, omega(0);
  Vector     p(resource(x)), phat(resource(x)), s(resource(x)), shat(resource(x)), 
             t(resource(x)), v(resource(x)), r(resource(x)), rtilde(resource(x));

  r = b - A * x;
  rtilde = r;

  while (! iter.finished(r)) {
    ++iter;
    rho_1 = dot(rtilde, r);
    MTL_THROW_IF(rho_1 == 0.0, unexpected_orthogonality());

    if (iter.first())
      p = r;
    else {
      MTL_THROW_IF(omega == 0.0, unexpected_orthogonality());
      beta = (rho_1 / rho_2) * (alpha / omega);
      p = r + beta * (p - omega * v);
    }
    phat = solve(M, p);
    v = A * phat;

    gamma = dot(rtilde, v);
    MTL_THROW_IF(gamma == 0.0, unexpected_orthogonality());

    alpha = rho_1 / gamma;
    s = r - alpha * v;
    
    if (iter.finished(s)) {
      x += alpha * phat;
      break;
    }
    shat = solve(M, s);
    t = A * shat;
    omega = dot(t, s) / dot(t, t);

    x += omega * shat + alpha * phat;
    r = s - omega * t;

    rho_2 = rho_1;    
  }
  return iter;
}

/// Solver class for BiCGStab method; right preconditioner ignored (prints warning if not identity)
/** Methods inherited from \ref base_solver. **/
template < typename LinearOperator, typename Preconditioner= pc::identity<LinearOperator>, 
	   typename RightPreconditioner= pc::identity<LinearOperator> >
class bicgstab_solver
  : public base_solver< bicgstab_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator >
{
    typedef base_solver< bicgstab_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator > base;
  public:
    /// Construct solver from a linear operator; generate (left) preconditioner from it
    explicit bicgstab_solver(const LinearOperator& A) : base(A), L(A) 
    {
	if (!pc::static_is_identity<RightPreconditioner>::value)
	    std::cerr << "Right Preconditioner ignored!" << std::endl;
    }

    /// Construct solver from a linear operator and (left) preconditioner
    bicgstab_solver(const LinearOperator& A, const Preconditioner& L) : base(A), L(L) 
    {
	if (!pc::static_is_identity<RightPreconditioner>::value)
	    std::cerr << "Right Preconditioner ignored!" << std::endl;
    }

    /// Solve linear system approximately as specified by \p iter
    template < typename HilbertSpaceX, typename HilbertSpaceB, typename Iteration >
    int solve(HilbertSpaceX& x, const HilbertSpaceB& b, Iteration& iter) const
    {
	return bicgstab(this->A, x, b, L, iter);
    }

  private:
    Preconditioner        L;
};

} // namespace itl

#endif // ITL_BICGSTAB_INCLUDE
