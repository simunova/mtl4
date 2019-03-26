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

// This test is written by Jan Bos to test the convergence of complex linear systems

#include <cmath>

// #include <boost/test/minimal.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <boost/numeric/itl/krylov/bicgstab_ell.hpp>


struct Random
{
  std::complex<double> operator()() const
  {
    return std::complex<double>(
        (2.*(static_cast<double>(rand())/RAND_MAX - 0.5)),
        (2.*(static_cast<double>(rand())/RAND_MAX - 0.5))
        );
  }
};

template <typename Matrix>
void test1(Matrix& m, double tau)
{
  Random m_rand;
  mtl::mat::inserter<Matrix> ins(m);
  size_t nrows=num_rows(m);
  std::complex<double> val;
  for (size_t r=0;r<nrows;++r)
  {
    for (size_t c=0;c<nrows;++c)
    {
      if(r==c) 
        ins(r,c) << 1.;
      else
      {
        val=m_rand();
        if (abs(val)<tau)
          ins(r,c) << val;
      }
    }
  }
}

int main(int, char**)
{
  const int N = 10; // Original from Jan had 2000 
  const int Niter = 3*N;

  using mtl::compressed2D; using mtl::dense_vector;
  using itl::cyclic_iteration; using itl::pc::identity; using itl::pc::ilu_0;

  typedef compressed2D<std::complex<double> > matrix_type;
  //typedef compressed2D<std::complex<double> ,mat::parameters<tag::col_major> > matrix_type;
  matrix_type                   A(N, N);
  dense_vector<std::complex<double> > b(N, std::complex<double>(1,0)), x(N);

  test1(A,0.194);
  std::cout << "A has " << A.nnz() << " non-zero entries" << std::endl;
  identity<matrix_type>     Ident(A);

  std::cout << "Non- preconditioned bicgstab" << std::endl;
  std::cout << "Won't convergence (for large examples)!" << std::endl;
  x= 0.5;
  cyclic_iteration<double> iter_0(b, Niter, 1.e-8, 0.0, 5);
  bicgstab(A, x, b, Ident, iter_0);

  std::cout << "Non-preconditioned bicgstab(1)" << std::endl;
  std::cout << "Won't convergence (for large examples)!" << std::endl;
  x= 0.5;
  cyclic_iteration<double> iter_1(b, Niter, 1.e-8, 0.0, 5);
  bicgstab_ell(A, x, b, Ident, Ident, iter_1,1);
 
  std::cout << "Non-preconditioned bicgstab(2)" << std::endl;
  x= 0.5;
  cyclic_iteration<double> iter_2b(b, Niter, 1.e-8, 0.0, 5);
  bicgstab_ell(A, x, b, Ident, Ident, iter_2b,2);

  std::cout << "Non-preconditioned bicgstab(4)" << std::endl;
  x= 0.5;
  cyclic_iteration<double> iter_4b(b, Niter, 1.e-8, 0.0, 5);
  bicgstab_ell(A, x, b, Ident, Ident, iter_4b,4);

  std::cout << "Non-preconditioned bicgstab(8)" << std::endl;
  x= 0.5;
  cyclic_iteration<double> iter_8b(b, Niter, 1.e-8, 0.0, 5);
  bicgstab_ell(A, x, b, Ident, Ident, iter_8b,8);

  ilu_0<matrix_type>        P(A);
  
  std::cout << "Right ilu(0) preconditioned bicgstab(1)" << std::endl;
  x= 0.5;
  cyclic_iteration<double> iter_1r(b, Niter, 1.e-8, 0.0, 5);
  bicgstab_ell(A, x, b, Ident, P, iter_1r,1);
 
  std::cout << "Right ilu(0) preconditioned bicgstab(2)" << std::endl;
  x= 0.5;
  cyclic_iteration<double> iter_2r(b, Niter, 1.e-8, 0.0, 5);
  bicgstab_ell(A, x, b, Ident, P, iter_2r,2);

  std::cout << "Left ilu(0) preconditioned bicgstab(4)" << std::endl;
  x= 0.5;
  cyclic_iteration<double> iter_4l(b, Niter, 1.e-8, 0.0, 5);
  bicgstab_ell(A, x, b, P, Ident, iter_4l,4);

  std::cout << "Right ilu(0) preconditioned bicgstab(4)" << std::endl;
  x= 0.5;
  cyclic_iteration<double> iter_4r(b, Niter, 1.e-8, 0.0, 5);
  bicgstab_ell(A, x, b, Ident, P, iter_4r,4);

  std::cout << "Right ilu(0) preconditioned bicgstab(8)" << std::endl;
  x= 0.5;
  cyclic_iteration<double> iter_8r(b, Niter, 1.e-8, 0.0, 5);
  bicgstab_ell(A, x, b, Ident, P, iter_8r,8);


  return 0;
}
