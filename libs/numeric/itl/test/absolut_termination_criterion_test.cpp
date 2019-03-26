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

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>


int main()
{
  // For a more realistic example set size to 1000 or larger
  const int size = 4, N = size * size;
  const double abs= 1.e-6; // Choose appropriately for size

  typedef mtl::compressed2D<double>  matrix_type;
  mtl::compressed2D<double>          A(N, N);
  laplacian_setup(A, size, size);

  itl::pc::diagonal<matrix_type>     P(A);

  mtl::dense_vector<double> x(N, 1.0), b(N);

  b = A * x;
  x= 0;

  itl::cyclic_iteration<double> iter(b, 500, 0, abs, 5);
  cg(A, x, b, P, iter);

  mtl::dense_vector<double> r(b - A * x);
  std::cout << "Norm r = " << two_norm(r) << "\n";
  if (two_norm(r) > abs) throw "Residuum not reduced enough.";

  return 0;
}
