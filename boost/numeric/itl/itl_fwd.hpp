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

#ifndef ITL_ITL_FWD_INCLUDE
#define ITL_ITL_FWD_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>

namespace itl {

    template <class Real>                class basic_iteration;
    template <class Real, class OStream> class cyclic_iteration;
    template <class Real, class OStream> class noisy_iteration;

    template <typename Solver, typename VectorIn, bool trans> class solver_proxy;

    namespace pc {

	template <typename PC, typename Vector, bool> struct solver;
	template <typename Matrix, typename Value> class identity;
	// template <typename Matrix, typename Value, typename Vector> Vector solve(const identity<Matrix>&, const Vector& x);	
	// template <typename Matrix, typename Vector> Vector adjoint_solve(const identity<Matrix>&, const Vector& x);

	template <typename Matrix, typename Value> class diagonal;
	// template <typename Matrix, typename Vector> Vector solve(const diagonal<Matrix>& P, const Vector& x);
	// template <typename Matrix, typename Vector> Vector adjoint_solve(const diagonal<Matrix>& P, const Vector& x);

	template <typename Matrix, typename Factorizer, typename Value> class ilu;
	template <typename Matrix, typename Value> class ilu_0; // Maybe we should declare the default here???
	template <typename Matrix, typename Value> class ilut; // Maybe we should declare the default here???
	// template <typename Matrix, typename Vector> Vector solve(const ilu_0<Matrix>& P, const Vector& x);
	// template <typename Matrix, typename Vector> Vector adjoint_solve(const ilu_0<Matrix>& P, const Vector& x);

	template <typename Matrix, typename Value> class ic_0; // Maybe we should declare the default here???
	// template <typename Matrix, typename Vector> Vector solve(const ic_0<Matrix>& P, const Vector& x);
	// template <typename Matrix, typename Value, typename Vector> Vector adjoint_solve(const ic_0<Matrix, Value>& P, const Vector& x);

    } //  namespace pc

    template < typename LinearOperator, typename HilbertSpaceX, typename HilbertSpaceB, 
	       typename Preconditioner, typename Iteration >
    int cg(const LinearOperator& A, HilbertSpaceX& x, const HilbertSpaceB& b, 
	   const Preconditioner& M, Iteration& iter);

    template < typename LinearOperator, typename Preconditioner= pc::identity<LinearOperator, double>, 
	       typename RightPreconditioner= pc::identity<LinearOperator, double> >
    class cg_solver;


    template < typename LinearOperator, typename Vector, 
	       typename Preconditioner, typename Iteration >
    int bicg(const LinearOperator &A, Vector &x, const Vector &b,
	     const Preconditioner &M, Iteration& iter);

    template < class LinearOperator, class HilbertSpaceX, class HilbertSpaceB, 
	       class Preconditioner, class Iteration >
    int bicgstab(const LinearOperator& A, HilbertSpaceX& x, const HilbertSpaceB& b, 
		 const Preconditioner& M, Iteration& iter);

    template < class LinearOperator, class HilbertSpaceX, class HilbertSpaceB, 
	       class Preconditioner, class Iteration >
    int bicgstab_2(const LinearOperator& A, HilbertSpaceX& x, const HilbertSpaceB& b, 
		   const Preconditioner& M, Iteration& iter);

    template < typename LinearOperator, typename Vector, 
	       typename LeftPreconditioner, typename RightPreconditioner, 
	       typename Iteration >
    int bicgstab_ell(const LinearOperator &A, Vector &x, const Vector &b,
		     const LeftPreconditioner &L, const RightPreconditioner &R, 
		     Iteration& iter, size_t l);

    template < typename LinearOperator, typename Preconditioner= pc::identity<LinearOperator, double>, 
	       typename RightPreconditioner= pc::identity<LinearOperator, double> >
    class bicgstab_ell_solver;

    template <typename Solver, unsigned N, bool Stored= false>
    class repeating_solver;

} // namespace itl

#endif // ITL_ITL_FWD_INCLUDE
