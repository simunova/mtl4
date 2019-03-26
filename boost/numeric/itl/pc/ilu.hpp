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

#ifndef ITL_PC_ILU_INCLUDE
#define ITL_PC_ILU_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/operation/adjoint.hpp>
#include <boost/numeric/mtl/operation/lower_trisolve.hpp>
#include <boost/numeric/mtl/operation/upper_trisolve.hpp>
#include <boost/numeric/mtl/operation/lu.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/itl/pc/solver.hpp>

namespace itl { namespace pc {

/// Incomplete LU factorization of a \p Matrix into matrices of type \p Value using a \p Factorizer (e.g. ILU(0) or ILUT)
template <typename Matrix, typename Factorizer, typename Value= typename mtl::Collection<Matrix>::value_type>
class ilu
{
  public:
    typedef Value                                         value_type;
    typedef typename mtl::Collection<Matrix>::size_type   size_type;
    typedef ilu                                           self;
    typedef Factorizer                                    factorizer_type;

    typedef mtl::mat::parameters<mtl::row_major, mtl::index::c_index, mtl::non_fixed::dimensions, false, size_type> para;
    typedef mtl::mat::compressed2D<value_type, para>                     L_type;
    typedef mtl::mat::compressed2D<value_type, para>                     U_type;
    typedef typename mtl::mat::traits::adjoint<L_type>::type     adjoint_L_type;
    typedef typename mtl::mat::traits::adjoint<U_type>::type     adjoint_U_type;

    typedef mtl::mat::detail::lower_trisolve_t<L_type, mtl::tag::unit_diagonal, true>    lower_solver_t;
    typedef mtl::mat::detail::upper_trisolve_t<U_type, mtl::tag::inverse_diagonal, true> upper_solver_t;

    typedef mtl::mat::detail::lower_trisolve_t<adjoint_U_type, mtl::tag::inverse_diagonal, true> adjoint_lower_solver_t;
    typedef mtl::mat::detail::upper_trisolve_t<adjoint_L_type, mtl::tag::unit_diagonal, true>    adjoint_upper_solver_t;

    /// Factorization adapted from Saad
    explicit ilu(const Matrix& A) 
      : f(A, L, U), lower_solver(L), upper_solver(U), adjoint_L(adjoint(L)), adjoint_U(adjoint(U)), 
	adjoint_lower_solver(adjoint_U), adjoint_upper_solver(adjoint_L) {}

    template <typename FactPara>
    ilu(const Matrix& A, const FactPara& p) 
      : f(A, p, L, U), lower_solver(L), upper_solver(U), adjoint_L(adjoint(L)), adjoint_U(adjoint(U)), 
	adjoint_lower_solver(adjoint_U), adjoint_upper_solver(adjoint_L) {}

    /// Solve  LU y = x --> y= U^{-1} L^{-1} x
    template <typename Vector>
    Vector solve(const Vector& x) const
    {
	Vector y;
	solve(x, y);
	return y;
    }

    // solve x = L y --> y0= L^{-1} x
    template <typename VectorIn, typename VectorOut>
    const VectorOut& solve_lower(const VectorIn& x, VectorOut&) const
    {
	static VectorOut y0;
	y0.change_resource(resource(x));
	lower_solver(x, y0);
	return y0;
    }

    // Solve  LU y = x --> y= U^{-1} L^{-1} x
    template <typename VectorIn, typename VectorOut>
    void solve(const VectorIn& x, VectorOut& y) const
    {
	mtl::vampir_trace<5039> tracer;
	const VectorOut& y0= solve_lower(x, y);

	y.checked_change_resource(x);
	upper_solver(y0, y);
    }


    /// Solve (LU)^H y = x --> y= L^{-H} U^{-H} x
    template <typename Vector>
    Vector adjoint_solve(const Vector& x) const
    {
	mtl::vampir_trace<5040> tracer;
	Vector y(resource(x));
	adjoint_solve(x, y);
	return y;
    }

    /// Solve (LU)^H y = x --> y= L^{-H} U^{-H} x
    template <typename VectorIn, typename VectorOut>
    void adjoint_solve(const VectorIn& x, VectorOut& y) const
    {
	mtl::vampir_trace<5040> tracer;
	y.checked_change_resource(x);
	// y= unit_upper_trisolve(adjoint(L), inverse_lower_trisolve(adjoint(U), x));
	static VectorOut y0;
	y0.change_resource(resource(x));
	adjoint_lower_solver(x, y0);
	adjoint_upper_solver(y0, y);
    }


    L_type get_L() { return L; }
    U_type get_U() { return U; }

  public:
    L_type                      L;
    U_type                      U;
  private:
    Factorizer                  f;
    lower_solver_t              lower_solver;
    upper_solver_t              upper_solver;
    adjoint_L_type              adjoint_L;
    adjoint_U_type              adjoint_U;
    adjoint_lower_solver_t      adjoint_lower_solver;
    adjoint_upper_solver_t      adjoint_upper_solver;

}; 

template <typename Value, typename Factorizer, typename V2>
class ilu<mtl::mat::dense2D<Value, mtl::mat::parameters<> >, Factorizer, V2> // last 2 arguments are dummies
{
  public:
    typedef mtl::mat::dense2D<Value, mtl::mat::parameters<> >    Matrix;
    typedef typename mtl::Collection<Matrix>::value_type  value_type;
    typedef typename mtl::Collection<Matrix>::size_type   size_type;
    typedef ilu                                           self;
    typedef Matrix                                        LU_type;

    ilu(const Matrix& A) : LU(A) { lu(LU, P); }

    // Solve  P^{-1}LU x = b --> x= U^{-1} L^{-1} P b
    template <typename Vector>
    Vector solve(const Vector& b) const { return lu_apply(LU, P, b); }

    // Solve  P^{-1}LU x = b --> x= U^{-1} L^{-1} P b
    template <typename VectorIn, typename VectorOut>
    void solve(const VectorIn& b, VectorOut& x) const { x= lu_apply(LU, P, b); }

    // Solve (P^{-1}LU)^H x = b --> x= P^{-1}L^{-H} U^{-H} b // P^{-1}^{-1}^H = P^{-1})
    template <typename Vector>
    Vector adjoint_solve(const Vector& b) const { return lu_adjoint_apply(LU, P, b); }

    // Solve (P^{-1}LU)^H x = b --> x= P^{-1}L^{-H} U^{-H} b // P^{-1}^{-1}^H = P^{-1})
    template <typename VectorIn, typename VectorOut>
    void adjoint_solve(const VectorIn& b, VectorOut& x) const { x= lu_adjoint_apply(LU, P, b); }

  private:
    LU_type                        LU;
    mtl::vec::dense_vector<size_type, mtl::vec::parameters<> >   P;
};

#if 0
template <typename Matrix, typename Factorizer, typename Value, typename Vector>
struct ilu_solver
  : mtl::assigner<ilu_solver<Matrix, Factorizer, Value, Vector> >
{
    typedef ilu<Matrix, Factorizer, Value> pc_type;

    ilu_solver(const pc_type& P, const Vector& x) : P(P), x(x) {}

    template <typename VectorOut>
    void assign_to(VectorOut& y) const
    {	P.solve(x, y);    }    

    const pc_type&       P; 
    const Vector&        x;
};
#endif



/// Solve LU x = b --> x= U^{-1} L^{-1} b
template <typename Matrix, typename Factorizer, typename Value, typename Vector>
// ilu_solver<Matrix, Factorizer, Value, Vector> 
solver<ilu<Matrix, Factorizer, Value>, Vector, false>
inline solve(const ilu<Matrix, Factorizer, Value>& P, const Vector& x)
{
    return solver<ilu<Matrix, Factorizer, Value>, Vector, false>(P, x);
}


/// Solve (LU)^H x = b --> x= L^{-H} U^{-H} b
template <typename Matrix, typename Factorizer, typename Value, typename Vector>
// Vector 
solver<ilu<Matrix, Factorizer, Value>, Vector, true>
inline adjoint_solve(const ilu<Matrix, Factorizer, Value>& P, const Vector& b)
{
    return solver<ilu<Matrix, Factorizer, Value>, Vector, true>(P, b);
}

// ic_0_evaluator not needed IC(0) and ILU(0) do the same at the upper triangle ;-)

}} // namespace itl::pc


#endif // ITL_PC_ILU_INCLUDE
