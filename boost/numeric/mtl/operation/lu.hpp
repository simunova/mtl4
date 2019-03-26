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

// With contributions from Cornelius Steinhardt

#ifndef MTL_MATRIX_LU_INCLUDE
#define MTL_MATRIX_LU_INCLUDE

#include <cmath>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/utility/lu_matrix_type.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/concept/magnitude.hpp>
#include <boost/numeric/mtl/matrix/upper.hpp>
#include <boost/numeric/mtl/matrix/lower.hpp>
#include <boost/numeric/mtl/matrix/permutation.hpp>
#include <boost/numeric/mtl/operation/adjoint.hpp>
#include <boost/numeric/mtl/operation/lower_trisolve.hpp>
#include <boost/numeric/mtl/operation/upper_trisolve.hpp>
#include <boost/numeric/mtl/operation/max_pos.hpp>
#include <boost/numeric/mtl/operation/permute.hpp>
#include <boost/numeric/mtl/operation/swap_row.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace mat {

/// LU factorization in place (without pivoting and optimization so far)
/** eps is tolerance for pivot element. If less or equal the matrix is considered singular.
    eps is given as double right now, might be refactored to the magnitude type of the value type in the future. **/
template <typename Matrix>
void inline lu(Matrix& LU, typename Magnitude<typename Collection<Matrix>::value_type>::type eps= 0)
{
    vampir_trace<5023> tracer;
    using std::abs;
    MTL_THROW_IF(num_rows(LU) != num_cols(LU), matrix_not_square());

    for (std::size_t k= 0; k < num_rows(LU)-1; k++) {
	if (abs(LU[k][k]) <= eps) throw matrix_singular(); 
	irange r(k+1, imax); // Interval [k+1, n-1]
	LU[r][k]/= LU[k][k];
	LU[r][r]-= LU[r][k] * LU[k][r];
    }
}

/// LU factorization in place (with pivoting and without optimization so far)
/** eps is tolerance for pivot element. If less or equal the matrix is considered singular.
    eps is given as double right now, might be refactored to the magnitude type of the value type in the future. **/
template <typename Matrix, typename PermuationVector>
typename mtl::traits::enable_if_vector<PermuationVector>::type
lu(Matrix& A, PermuationVector& P, typename Magnitude<typename Collection<Matrix>::value_type>::type eps= 0)
{
    vampir_trace<5024> tracer;
    using math::zero; using std::abs;
    typedef typename Collection<Matrix>::size_type    size_type;
	typedef typename Collection<PermuationVector>::value_type value_p_type;
    size_type ncols = num_cols(A), nrows = num_rows(A);

    MTL_THROW_IF(ncols != nrows , matrix_not_square());
    P.change_dim(nrows);

    for (value_p_type i= 0; i < value_p_type(nrows); i++)
        P[i]= i;

    for (size_type i= 0; i < nrows-1; i++) {
	irange r(i+1, imax), ir(i, i+1); // Intervals [i+1, n-1], [i, i]
	size_type rmax= max_abs_pos(A[irange(i, imax)][ir]).first + i;
	swap_row(A, i, rmax); 
	swap_row(P, i, rmax);
	
	if(abs(A[i][i]) <= eps) throw matrix_singular(); // other gmres test doesn't work
       
	A[r][i]/= A[i][i];              // Scale column i
	A[r][r]-= A[r][i] * A[i][r]; 	 // Decrease bottom right block of matrix
    }
}


/// LU factorization without factorization that returns the matrix
/** eps is tolerance for pivot element. If less or equal the matrix is considered singular.
    eps is given as double right now, might be refactored to the magnitude type of the value type in the future. **/
template <typename Matrix>
Matrix inline lu_f(const Matrix& A, typename Magnitude<typename Collection<Matrix>::value_type>::type eps= 0)
{
    vampir_trace<5025> tracer;
    Matrix LU(A);
    lu(LU, eps);
    return LU;
}

/// Solve Ax = b by LU factorization without pivoting; vector x is returned
template <typename Matrix, typename Vector>
Vector inline lu_solve_straight(const Matrix& A, const Vector& b, typename Magnitude<typename Collection<Matrix>::value_type>::type eps= 0)
{
    vampir_trace<5026> tracer;
    Matrix LU(A);
    lu(LU, eps);
    return upper_trisolve(upper(LU), unit_lower_trisolve(strict_lower(LU), b));
}

/// Solve LUx = b by with forward and backward-LU subsitution, lu(LU) was allready done
template <typename Matrix, typename Vector>
Vector inline lu_solve_apply(const Matrix& LU, const Vector& b)
{
    vampir_trace<5026> tracer;
    return upper_trisolve(upper(LU), unit_lower_trisolve(strict_lower(LU), b));
}

/// Apply the factorization L*U with permutation P on vector b to solve Ax = b
template <typename Matrix, typename PermVector, typename Vector>
Vector inline lu_apply(const Matrix& LU, const PermVector& P, const Vector& b)
{
    vampir_trace<5027> tracer;
    return upper_trisolve(upper(LU), unit_lower_trisolve(strict_lower(LU), Vector(permute(P, b))));
}


/// Solve Ax = b by LU factorization with column pivoting; vector x is returned
template <typename Matrix, typename Vector>
Vector inline lu_solve(const Matrix& A, const Vector& b, typename Magnitude<typename Collection<Matrix>::value_type>::type eps= 0)
{
    vampir_trace<5028> tracer;
    mtl::dense_vector<std::size_t, vec::parameters<> > P(num_rows(A));
    Matrix                    LU(A);

    lu(LU, P, eps);
    return lu_apply(LU, P, b);
}


/// Apply the factorization L*U with permutation P on vector b to solve adjoint(A)x = b
/** That is \f$P^{-1}(LU)^H x = b\f$ --> \f$x= P^{-1}L^{-H} U^{-H} b\f$ where \f$P^{{-1}^{{-1}^H}} = P^{-1}\f$ **/
template <typename Matrix, typename PermVector, typename Vector>
Vector inline lu_adjoint_apply(const Matrix& LU, const PermVector& P, const Vector& b)
{
    vampir_trace<5029> tracer;
    return Vector(reverse_permute(P, unit_upper_trisolve(adjoint(LU), lower_trisolve(adjoint(LU), b))));
}


/// Solve \f$adjoint(A)x = b\f$ by LU factorization with column pivoting; vector x is returned
template <typename Matrix, typename Vector>
Vector inline lu_adjoint_solve(const Matrix& A, const Vector& b, typename Magnitude<typename Collection<Matrix>::value_type>::type eps= 0)
{
    vampir_trace<5030> tracer;
    mtl::dense_vector<std::size_t, vec::parameters<> > P(num_rows(A));
    Matrix                    LU(A);

    lu(LU, P, eps);
    return lu_adjoint_apply(LU, P, b);
}

/// Class that keeps LU factorization (and permutation); using column pivoting
template <typename Matrix>
class lu_solver
{
    typedef typename mtl::traits::lu_matrix_type<Matrix>::type            matrix_type;
    typedef mtl::vec::dense_vector<std::size_t, mtl::vec::parameters<> > permutation_type;
  public:
    /// Construct from matrix \p A and use optionally threshold \p eps in factorization
    explicit lu_solver(const Matrix& A, typename Magnitude<typename Collection<Matrix>::value_type>::type eps= 0) 
      : LU(A), P(num_rows(A))
    {
	lu(LU, P, eps);
    }

    /// Solve A*x = b with factorization from constructor
    template <typename VectorIn, typename VectorOut>
    void solve(const VectorIn& b, VectorOut& x) const
    {
	x= upper_trisolve(upper(LU), unit_lower_trisolve(strict_lower(LU), VectorIn(permute(P, b))));
    }
    /// Solve \f$adjoint(A)x = b\f$ using LU factorization
    template <typename VectorIn, typename VectorOut>
    void adjoint_solve(const VectorIn& b, VectorOut& x) const
    {
	x= reverse_permute(P, unit_upper_trisolve(adjoint(LU), lower_trisolve(adjoint(LU), b)));
    }

  private:
    matrix_type      LU;
    permutation_type P;
};


}} // namespace mtl::matrix 

















// ### For illustration purposes
#if 0

namespace mtl { namespace matrix {

template <typename Matrix>
void inline lu(Matrix& LU)
{
    MTL_THROW_IF(num_rows(LU) != num_cols(LU), matrix_not_square());

    typedef typename Collection<Matrix>::value_type   value_type;
    typedef typename Collection<Matrix>::size_type    size_type;

    size_type n= num_rows(LU);
    for (size_type k= 0; k < num_rows(LU); k++) {
	value_type pivot= LU[k][k];
	for (size_type j= k+1; j < n; j++) {
	    value_type alpha= LU[j][k]/= pivot;
	    for (size_type i= k+1; i < n; i++)
		LU[j][i]-= alpha * LU[k][i];
	}
    }
}



}} // namespace mtl::matrix

#endif

#endif // MTL_MATRIX_LU_INCLUDE
