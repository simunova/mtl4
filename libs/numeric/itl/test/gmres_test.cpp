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


int main(int, char**)
{



  const int N = 2; // Original from Jan had 2000
  const int Niter = 5000;
  using itl::pc::identity; using itl::pc::ilu_0; using itl::pc::ic_0; using itl::pc::diagonal;
  //typedef mtl::dense2D<double> matrix_type;
  //typedef compressed2D<std::complex<double>mat::parameters<tag::col_major> > matrix_type;
  typedef mtl::compressed2D<double> matrix_type;
  matrix_type                   A(N, N);
  laplacian_setup(A, N, N);
  mtl::dense_vector<double> b(N*N, 1), x(N*N),r(N*N);
  identity<matrix_type>         Ident(A);
  //ic_0<matrix_type>             ic(A);
  //ilu_0<matrix_type>            ilu(A);
  //diagonal<matrix_type>         diag(A);

  //test1(A,0.194);
  std::cout << "A has " << A.nnz() << " non-zero entries" << std::endl;
  std::cout << "A =\n"  << A << " \n";

  std::cout << "Non- preconditioned bicgstab  Won't convergence (for large examples)!" << std::endl;
  x= 2.0;
  itl::basic_iteration<double> iter_0(b, Niter, 1.e-8);

  bicgstab(A, x, b, Ident, iter_0);
  r= A*x-b;

  std::cout << "START GMRES  Won't convergence (for large examples)!" << std::endl;
  std::cout << "\n Non-preconditioned gmres(1)" << std::endl;
  x= 5.0;
  itl::basic_iteration<double> iter_1(b, 1, 1.e-8);
  gmres(A, x, b, Ident, Ident, iter_1, 1);
  r= A*x-b;
  if (two_norm(r) > 0.00005) throw "gmres doesn't converge";

  std::cout << "\n Non-preconditioned gmres(2) (doesn't converge even for the test)" << std::endl;
  x= 2.0;
  itl::basic_iteration<double> iter_2(b, 8, 1.e-8);
  gmres(A, x, b, Ident, Ident, iter_2, 2);
  r= A*x-b;
  // if (two_norm(r) > 0.00005) throw "gmres(2) doesn't converge";
  if (two_norm(r) > 0.00005) std::cout << "GMRES(2) didn't converge after 8 titerations.\n";

#if 1
  std::cout << "\n Non-preconditioned gmres(4)" << std::endl;
  x= 2.5;
  itl::basic_iteration<double> iter_4(b, 16, 1.e-8);
  gmres(A, x, b, Ident, Ident,  iter_4, 4);
  r= A*x-b;
  // if (two_norm(r) > 0.000001) throw "gmres(4) doesn't converge even with more iterations and restarts";
  if (two_norm(r) > 0.00005) std::cout << "GMRES(4) didn't converge after 16 titerations.\n";

  std::cout << "\n Non-preconditioned gmres(4) more iterations " << std::endl;
  x= 2.5;
  itl::basic_iteration<double> iter_5(b, 32, 1.e-8);
  gmres(A, x, b, Ident, Ident,  iter_5, 4);
  r= A*x-b;
  // if (two_norm(r) > 0.00005) throw "gmres(4) doesn't converge";

  std::cout << "\n Non-preconditioned gmres(8)" << std::endl;
  x= 2.5;
  itl::basic_iteration<double> iter_8(b, 8, 1.e-8);
  gmres(A, x, b, Ident, Ident, iter_8, 8);
  r= A*x-b;
  if (two_norm(r) > 0.00005) throw "gmres(8) doesn't converge";

  std::cout << "\n Non-preconditioned gmres(16)" << std::endl;
  x= 2.5;
  itl::basic_iteration<double> iter_16(b, 32, 1.e-8);
  gmres(A, x, b, Ident, Ident, iter_16, 16);
  r= A*x-b;
  if (two_norm(r) > 0.00005) throw "gmres(16) doesn't converge";

  std::cout << "\n Non-preconditioned gmres(32)" << std::endl;
  x= 2.5;
  itl::basic_iteration<double> iter_32(b, 32, 1.e-8);
  gmres(A, x, b, Ident, Ident, iter_32, 32);
  r= A*x-b;
  if (two_norm(r) > 0.00005) throw "gmres(32) doesn't converge";
#endif

  test1(A,0.194);
  std::cout << "A has " << A.nnz() << " non-zero entries" << std::endl;
  std::cout << "A =\n"  << A << " \n";

  return 0;
}
