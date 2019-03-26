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

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

template <typename Matrix>
void fill(Matrix& m)
{
  double s=19, u=21, p=16, e=5, r=18, l =12;

  std::complex<double> delta(2,1);
  mtl::mat::inserter<Matrix> sm(m);
  // set diagonal
  sm(0,0) << s;
  sm(1,1) << u;
  sm(2,2) << p;
  sm(3,3) << e;
  sm(4,4) << r;
  // below diagonal
  sm(1,0) << l+delta;
  sm(4,0) << l;
  sm(2,1) << l;
  sm(4,1) << l+delta;
  // above diagonal
  sm(0,2) << u;
  sm(0,3) << u+delta;
  //sm(3,4) << u;
  sm(3,4) << u + delta;
}

int main()
{
  // For a more realistic example set size to 1000 or larger
  const int N = 5;

  typedef mtl::compressed2D<std::complex<double> > matrix_type;
  matrix_type                   A(N, N);

  itl::pc::identity<matrix_type>     P(A);
  mtl::dense_vector<std::complex<double> > b(N, std::complex<double>(1,1)), x(N);

  fill(A);
  x= 0;

  itl::cyclic_iteration<double> iter(b, N, 1.e-6, 0.0, 5);
  bicgstab(A, x, b, P, iter);
 
  return 0;
}
