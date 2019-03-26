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

#include <cmath>
#include <string>

// #include <boost/test/minimal.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>




template <typename Matrix>
void test1(Matrix& m, double tau)
{
  mtl::mat::inserter<Matrix> ins(m);
  size_t nrows=num_rows(m);
  double val;
  for (size_t r=0;r<nrows;++r)
  {
    for (size_t c=0;c<nrows;++c)
    {
      if(r==c)
        ins(r,c) << 1.;
      else
      {
        val=2.*(static_cast<double>(rand())/RAND_MAX - 0.5);
        if (val<tau)
          ins(r,c) << val;
      }
    }
  }
}

template <typename Matrix, typename Vector, typename Left, typename Right>
void test(char const* name, char const* comment, Matrix const& A, Vector& x, Vector const& b, Left const& L, Right const& R, 
	  unsigned restart, bool check_convergence= true)
{
    const int Niter = 100;
    
    std::cout << name << comment << "\n";
    x= 2.0, 3., 4., 8;

    itl::cyclic_iteration<double> iter(b, Niter, 1.e-8, 0.0, 10);
    gmres(A, x, b, L, R, iter, restart);
    std::cout << "x= " << x << " \n" ;
    Vector r(b - A*x);
    if (false && check_convergence && two_norm(r) > 0.00001) 
	throw std::string(name) + std::string(" doesn't converge!");
}


int main(int, char**)
{
    const int N = 2;
    typedef mtl::compressed2D<double> matrix_type;
    matrix_type                   A(N*N, N*N);
    laplacian_setup(A, N, N);

    mtl::dense_vector<double> b(N*N, 1), x(N*N,1), r(N*N);
 
    itl::pc::identity<matrix_type>         Ident(A);
    itl::pc::ic_0<matrix_type>             ic(A);
    itl::pc::ilu_0<matrix_type>            ilu(A);
    itl::pc::diagonal<matrix_type>         diag(A);

    std::cout << "A has " << A.nnz() << " non-zero entries" << std::endl;
    std::cout << "A =\n" << A << " \n";

    test("Non-preconditioned GMRES(1)", "\nWon't convergence (for large examples,without restarts)!",
	 A, x, b, Ident, Ident, 1, false);
    test("Non-preconditioned GMRES(4)", "", A, x, b, Ident, Ident, 4);
    test("Left ILU(0) GMRES(4)", "", A, x, b, ilu, Ident, 4);
    test("Left IC(0) GMRES(4)", "", A, x, b, ic, Ident, 4);
    test("Left diag GMRES(4)", "", A, x, b, diag, Ident, 4);

    test("Right ILU(0) GMRES(4)", "", A, x, b, Ident, ilu, 4);
    test("Right IC(0) GMRES(4)", "", A, x, b, Ident, ic, 4);
    test("Right diag GMRES(4)", "", A, x, b, Ident, diag, 4);

    test("Left ILU(0) Right ILU(0) GMRES(4)", "", A, x, b, ilu, ilu, 4);
    test("Left ILU(0) Right IC(0) GMRES(4)", "", A, x, b, ilu, ic, 4);
    test("Left ILU(0) Right diag GMRES(4)", "", A, x, b, ilu, diag, 4);

    return 0;
}

