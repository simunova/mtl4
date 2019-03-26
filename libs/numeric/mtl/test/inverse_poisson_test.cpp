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

#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
 
using namespace std;  

typedef mtl::compressed2D<double> m_type;

void assemble_poisson2D(m_type& A, int n)
{
    m_type dia_block(n, n), off_block(n, n);
    {
	mtl::mat::inserter<m_type> ins(dia_block);
	for (int i= 0; i < n; i++) {
	    ins[i][i] << 4.0;
	    if (i > 0) ins[i][i-1] << -1.0;
	    if (i < n-1) ins[i][i+1] << -1.0;
	}
    }
    off_block= -1.0;
    A.change_dim(n*n, n*n);
    {
	mtl::mat::inserter<m_type> ins(A);
	for (int i= 0; i < n; i++) {
	    ins[i*n][i*n] << dia_block;
	    if (i > 0) ins[i*n][i*n-n] << off_block;
	    if (i < n-1) ins[i*n][i*n+n] << off_block;
	}
    }
}

//
// Inverse iteration: is almost a literal copy from algorithm 9.1.2
//
// LinearSolver: BinaryFunction
// EigenVector: GLAS dense vector
template <class LinearSolver, class EigenVector>
double inverse_iteration( LinearSolver& solver, EigenVector& x, int m ) 
{
    EigenVector y( size(x), 0.0 ) ;
    double lambda = 0.0 ;
    x /= two_norm(x) ;
    for (int i=0; i<m; ++i) {
	solver( x, y ) ;
	lambda = mtl::sum(x) / mtl::sum(y) ;
	x = y / two_norm(y) ;
    }
    return lambda ;
}

//
// Inverse iteration: is almost a literal copy from algorithm 9.1.2
//
// LinearSolver: BinaryFunction
// EigenVector: GLAS dense vector
template <class LinearSolver, class Matrix, class EigenVector>
double condition( LinearSolver& solver, const Matrix& A, EigenVector& x, int m ) {
  EigenVector y( size(x) ) ;
  double lambda= 0.0 ;
  x/= two_norm(x);
  for (int i=0; i<m; ++i) {
      y= A * x; 
      lambda = mtl::sum(y) / mtl::sum(x) ;
      x = y / two_norm(y) ;
  }
  return std::abs(lambda / inverse_iteration(solver, x, m));
}

template <typename Matrix>
class cg_solver
{
public:
    cg_solver(const Matrix& A) : A(A) {}

    template <typename VectorIn, typename VectorOut>
    void operator() (const VectorIn& v_in, VectorOut& v_out)
    {
	itl::pc::diagonal<Matrix>     P(A);
	itl::cyclic_iteration<double>  iter(v_in, 500, 1.e-6);

	cg(A, v_out, v_in, P, iter);
    }

private:
    const Matrix& A;
};

int main(int, char**)
{
    m_type A;
    assemble_poisson2D(A, 10);

    mtl::dense_vector<double>   x(num_rows(A));
    mtl::seed<double>           seed;
    random(x, seed);

    cg_solver<m_type> sol(A);
    double lambda= inverse_iteration(sol, x, 10);
    cout << "lambda is " << lambda << '\n';

    double cond= condition(sol, A, x, 10);
    cout << "cond is " << cond << '\n';

    return 0;
}
