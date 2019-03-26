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


#ifndef ITL_GMRES_INCLUDE
#define ITL_GMRES_INCLUDE

#include <algorithm>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/multi_vector.hpp>
#include <boost/numeric/mtl/operation/givens.hpp>
#include <boost/numeric/mtl/operation/two_norm.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>

#include <boost/numeric/itl/krylov/base_solver.hpp>
#include <boost/numeric/itl/pc/identity.hpp>

namespace itl {

/// Generalized Minimal Residual method (without restart)
/** It computes at most kmax_in iterations (or size(x) depending on what is smaller) 
    regardless on whether the termination criterion is reached or not.   **/
template < typename Matrix, typename Vector, typename LeftPreconditioner, typename RightPreconditioner, typename Iteration >
int gmres_full(const Matrix &A, Vector &x, const Vector &b,
               LeftPreconditioner &L, RightPreconditioner &R, Iteration& iter)
{
    using mtl::size; using mtl::irange; using mtl::iall; using std::abs; using std::sqrt;
    typedef typename mtl::Collection<Vector>::value_type Scalar;
    typedef typename mtl::Collection<Vector>::size_type  Size;

    if (size(b) == 0) throw mtl::logic_error("empty rhs vector");

    const Scalar                zero= math::zero(Scalar());
    Scalar                      rho, nu, hr;
    Size                        k, kmax(std::min(size(x), Size(iter.max_iterations() - iter.iterations())));
    Vector                      r0(b - A *x), r(solve(L,r0)), va(resource(x)), va0(resource(x)), va00(resource(x));
    mtl::mat::multi_vector<Vector>   V(Vector(resource(x), zero), kmax+1); 
    mtl::dense_vector<Scalar>   s(kmax+1, zero), c(kmax+1, zero), g(kmax+1, zero), y(kmax, zero);  // replicated in distributed solvers 
    mtl::mat::dense2D<Scalar>        H(kmax+1, kmax);                                             // dito
    H= 0;

    rho= g[0]= two_norm(r);
    if (iter.finished(rho))
	return iter;
    V.vector(0)= r / rho;
    H= zero;

    // GMRES iteration
    for (k= 0; k < kmax ; ++k, ++iter) {
        va0= A * Vector(solve(R, V.vector(k)));
        V.vector(k+1)= va= solve(L,va0);
	// orth(V, V[k+1], false); 
        // modified Gram Schmidt method
        for (Size j= 0; j < k+1; j++) {
	    H[j][k]= dot(V.vector(j), V.vector(k+1));
	    V.vector(k+1)-= H[j][k] * V.vector(j);
        }
        H[k+1][k]= two_norm(V.vector(k+1));
        //reorthogonalize
        for(Size j= 0; j < k+1; j++) {
	    hr= dot(V.vector(k+1), V.vector(j));
            H[j][k]+= hr;
            V.vector(k+1)-= hr * V.vector(j);
        }
        H[k+1][k]= two_norm(V.vector(k+1));
	if (H[k+1][k] != zero)                // watch for breakdown    
            V.vector(k+1)*= 1. / H[k+1][k];

        // k Given's rotations
	for(Size i= 0; i < k; i++)
	    mtl::mat::givens<mtl::mat::dense2D<Scalar> >(H, H[i][k-1], H[i+1][k-1]).trafo(i);
	
       nu= sqrt(H[k][k]*H[k][k]+H[k+1][k]*H[k+1][k]);
       if(nu != zero){
            c[k]=  H[k][k]/nu;
            s[k]= -H[k+1][k]/nu;
            H[k][k]=c[k]*H[k][k]-s[k]*H[k+1][k];
            H[k+1][k]=0;
 	    mtl::vec::givens<mtl::vec::dense_vector<Scalar> >(g, c[k], s[k]).trafo(k);
        }
	rho= abs(g[k+1]);
    }
    
    //reduce k, to get regular matrix
    while (k > 0 && abs(g[k-1])<= iter.atol()) k--;

    // iteration is finished -> compute x: solve H*y=g as far as rank of H allows
    irange                  range(k);
    for (; !range.empty(); --range) {
	try {
	    y[range]= lu_solve(H[range][range], g[range]); 
	} catch (mtl::matrix_singular) { continue; } // if singular then try with sub-matrix
	break;
    }

    if (range.finish() < k)
  	std::cerr << "GMRES orhogonalized with " << k << " vectors but matrix singular, can only use " 
		  << range.finish() << " vectors!\n";
    if (range.empty())
        return iter.fail(2, "GMRES did not find any direction to correct x");
    x+= Vector(solve(R, Vector(V.vector(range)*y[range])));
    
    r= b - A*x;
    return iter.terminate(r);
}

/// Generalized Minimal Residual method with restart
template < typename Matrix, typename Vector, typename LeftPreconditioner,
           typename RightPreconditioner, typename Iteration >
int gmres(const Matrix &A, Vector &x, const Vector &b,
          LeftPreconditioner &L, RightPreconditioner &R,
	  Iteration& iter, typename mtl::Collection<Vector>::size_type restart)
{   
     do {
	 Iteration inner(iter);
	 inner.set_max_iterations(std::min(int(iter.iterations()+restart), iter.max_iterations()));
	 inner.suppress_resume(true);
	 gmres_full(A, x, b, L, R, inner);
	 iter.update_progress(inner);
     } while (!iter.finished());

     return iter;
}

/// Solver class for GMRES; right preconditioner ignored (prints warning if not identity)
/** Methods inherited from \ref base_solver. **/
template < typename LinearOperator, typename Preconditioner= pc::identity<LinearOperator>, 
	   typename RightPreconditioner= pc::identity<LinearOperator> >
class gmres_solver
  : public base_solver< gmres_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator >
{
    typedef base_solver< gmres_solver<LinearOperator, Preconditioner, RightPreconditioner>, LinearOperator > base;
  public:
    /// Construct solver from a linear operator; generate (left) preconditioner from it
    explicit gmres_solver(const LinearOperator& A, size_t restart= 8) 
      : base(A), restart(restart), L(A), R(A) {}

    /// Construct solver from a linear operator and left preconditioner
    gmres_solver(const LinearOperator& A, size_t restart, const Preconditioner& L) 
      : base(A), restart(restart), L(L), R(A) {}

    /// Construct solver from a linear operator and left preconditioner
    gmres_solver(const LinearOperator& A, size_t restart, const Preconditioner& L, const RightPreconditioner& R) 
      : base(A), restart(restart), L(L), R(R) {}

    /// Solve linear system approximately as specified by \p iter
    template < typename HilbertSpaceX, typename HilbertSpaceB, typename Iteration >
    int solve(HilbertSpaceX& x, const HilbertSpaceB& b, Iteration& iter) const
    {
	return gmres(this->A, x, b, L, R, iter, restart);
    }

  private:
    size_t                restart;
    Preconditioner        L;
    RightPreconditioner   R;
};


} // namespace itl

#endif // ITL_GMRES_INCLUDE


