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

// Peter Sonneveld and Martin B. van Gijzen, IDR(s): a family of simple and fast algorithms for solving large nonsymmetric linear systems. 
// SIAM J. Sci. Comput. Vol. 31, No. 2, pp. 1035-1062 (2008). (copyright SIAM)

#ifndef ITL_IDR_S_INCLUDE
#define ITL_IDR_S_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/operation/random.hpp>
#include <boost/numeric/mtl/operation/orth.hpp>
#include <boost/numeric/mtl/operation/resource.hpp>
#include <boost/numeric/mtl/matrix/strict_upper.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#include <boost/numeric/itl/krylov/base_solver.hpp>

namespace itl {

/// Induced Dimension Reduction on s dimensions (IDR(s)) 
template < typename LinearOperator, typename Vector, 
	   typename LeftPreconditioner, typename RightPreconditioner, 
	   typename Iteration >
int idr_s(const LinearOperator &A, Vector &x, const Vector &b,
	  const LeftPreconditioner &, const RightPreconditioner &, 
	  Iteration& iter, size_t s)
{
    mtl::vampir_trace<7010> tracer;
    using mtl::size; using mtl::iall; using mtl::mat::strict_upper;
    typedef typename mtl::Collection<Vector>::value_type Scalar;
    typedef typename mtl::Collection<Vector>::size_type  Size;

    if (size(b) == 0) throw mtl::logic_error("empty rhs vector");
    if (s < 1) s= 1;

    const Scalar                zero= math::zero(Scalar());
    Scalar                      omega(zero);
    Vector                      x0(x), y(resource(x)), v(resource(x)), t(resource(x)), q(resource(x)), r(b - A * x);
    mtl::mat::multi_vector<Vector>   dR(Vector(resource(x), zero), s), dX(Vector(resource(x), zero), s), P(Vector(resource(x), zero), s);
    mtl::dense_vector<Scalar>   m(s), c(s), dm(s);   // replicated in distributed solvers 
    mtl::mat::dense2D<Scalar>        M(s, s);             // dito

    random(P); 
    P.vector(0)= r;
    orth(P);

    for (size_t k= 0; k < s; k++) {
	v= A * r;
	omega= dot(v, r) / dot(v, v);
	dX.vector(k)= omega * r;
	dR.vector(k)= -omega * v;
	x+= dX.vector(k); 
	r+= dR.vector(k);
	if ((++iter).finished(r)) return iter;
	M[iall][k]= trans(P) * dR.vector(k); 
    }

    Size oldest= 0;
    m= trans(P) * r;

    while (! iter.finished(r)) {
       
	for (size_t k= 0; k < s; k++) {
	    c= lu_solve(M, m);
	    q= dR * -c;    
	    v= r + q;
	    if (k == 0) {
		t= A * v;
		omega= dot(t, v) / dot(t, t);
		dR.vector(oldest)= q - omega * t;
		dX.vector(oldest)= omega * v - dX * c;
	    } else {
		dX.vector(oldest)= omega * v - dX * c;
		dR.vector(oldest)= A * -dX.vector(oldest);
	    }
	    r+= dR.vector(oldest);
	    x+= dX.vector(oldest);

	    if ((++iter).finished(r))
		return iter;

	    dm= trans(P) * dR.vector(oldest);
	    M[iall][oldest]= dm;
	    m+= dm;
	    oldest= (oldest + 1) % s;
	}
    }
    return iter;
}

/// Solver class for IDR(s) method; right preconditioner ignored (prints warning if not identity)
/** Methods inherited from \ref base_solver. **/
template < typename LinearOperator, typename Preconditioner= pc::identity<LinearOperator>, 
	   typename RightPreconditioner= pc::identity<LinearOperator> >
class idr_s_solver
  : public base_solver< idr_s_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator >
{
    typedef base_solver< idr_s_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator > base;
  public:
  public:
    /// Construct solver from a linear operator; generate (left) preconditioner from it
    explicit idr_s_solver(const LinearOperator& A, size_t s= 8) : base(A), s(s), L(A), R(A) {}

    /// Construct solver from a linear operator and left preconditioner
    idr_s_solver(const LinearOperator& A, size_t s, const Preconditioner& L) : base(A), s(s), L(L), R(A) {}

    /// Construct solver from a linear operator and left preconditioner
    idr_s_solver(const LinearOperator& A, size_t s, const Preconditioner& L, const RightPreconditioner& R) 
      : base(A), s(s), L(L), R(R) {}

    /// Solve linear system approximately as specified by \p iter
    template < typename HilbertSpaceX, typename HilbertSpaceB, typename Iteration >
    int solve(HilbertSpaceX& x, const HilbertSpaceB& b, Iteration& iter) const
    {
	return idr_s(this->A, x, b, L, R, iter, s);
    }

  private:
    size_t                s;
    Preconditioner        L;
    RightPreconditioner   R;
};


} // namespace itl

#endif // ITL_IDR_S_INCLUDE
