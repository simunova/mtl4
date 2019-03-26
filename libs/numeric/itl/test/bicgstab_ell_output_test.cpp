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
#include <boost/numeric/mtl/operation/trans.hpp>

template <typename Matrix>
void test1(Matrix& m, double tau)
{
    mtl::mat::inserter<Matrix> ins(m);
    size_t nrows=num_rows(m);
    double val;
    for (size_t r= 0; r < nrows; ++r) 
	for (size_t c= 0; c < nrows; ++c) 
	    if (r == c)
		ins(r,c) << 1.;
	    else {
		val= 2.*(static_cast<double>(rand())/RAND_MAX - 0.5);
		if (val < tau)
		    ins(r,c) << val;
	    }			
}


int main(int, char**)
{

  const int N = 10; 
  const int Niter = 10*N;

  using itl::pc::identity; using itl::pc::ilu_0; using itl::pc::ic_0; using itl::pc::diagonal;
  //typedef mtl::dense2D<double> matrix_type;
  typedef mtl::compressed2D<double> matrix_type;
  matrix_type                   A(N, N);
  mtl::dense_vector<double>     b(N*N, 1), x(N*N), r(x);
  laplacian_setup(A, N, N);
  identity<matrix_type>         Ident(A);

  x= 0.5;
  itl::cyclic_iteration<double> iter_1(b, Niter, 0, 1.e-8, 10);

  bicgstab_ell(A, x, b, Ident, Ident, iter_1, 8);
  r= A*x-b;
  std::cout << "|A*x-b|=" << two_norm(r) << "\n";
  if (two_norm(r) > 0.000001) throw "bicgstab_ell doesn't converge";


  return 0;
}




