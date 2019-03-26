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

using namespace std;
int main()
{
  // For a more realistic example set size to 1000 or larger
  const int size = 8, N = size * size;

  //typedef mtl::dense2D<double>  matrix_type;
  typedef mtl::compressed2D<double>  matrix_type;
  matrix_type                      A(N, N), dr(5,5);

  laplacian_setup(A, size, size);

  dr= 1,  1,  1,   0,  0,
      1, -1, -2,   0,  0,
      1, -2,  1,   0,  0,
      0,  0,  0, -10,  0,
      0,  0,  0,   0, 22;

  itl::pc::diagonal<matrix_type>     P(A), Pb(dr);
  mtl::dense_vector<double>          x(N, 1.0), b(N), xb(5, 1.0), bb(5);
  mtl::dense_vector<complex<double> > xz(5,complex<double>(1.0, 0.0)), bz(5);

  bb= dr * xb;
  xb= 0;
  
  b= A * x;
  x= 0;

  itl::cyclic_iteration<double> iter(b, 500, 1.e-6, 0.0, 1);
  qmr(A, x, b, P, P, iter);
  std::cout<< "x=" << x << "\n";

  itl::cyclic_iteration<double> iterb(bb, 500, 1.e-6, 0.0, 1);
  qmr(dr, xb, bb, Pb, Pb, iterb);
  std::cout<< "xb=" << xb << "\n";

  return 0;
}


