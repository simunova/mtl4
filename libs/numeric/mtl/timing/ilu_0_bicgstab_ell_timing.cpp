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
// #define CRS_CVEC_MULT_NO_ACCEL

// #define CRS_CVEC_MULT_NO_ACCEL // for benchmarking
// #define MTL_LAZY_LOOP_WO_UNROLL // for benchmarking

#include <iostream>
#include <typeinfo>
#include <boost/timer.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>


int main()
{
    mtl::vampir_trace<9999> tracer;
  // For a more realistic example set size to 1000 or larger
  const int size = 1000, N = size * size; 
  using namespace mtl;

  typedef unsigned size_type;
  // typedef std::size_t size_type;
  std::cout << "sizeof in size_type is " << sizeof(size_type) << '\n';
  typedef mat::parameters<row_major, mtl::index::c_index, non_fixed::dimensions, false, size_type> para;
  typedef compressed2D<double, para>  matrix_type;
  matrix_type          A(N, N);
  laplacian_setup(A, size, size);

  itl::pc::ilu_0<matrix_type, float>     P(A);
  itl::pc::identity<matrix_type>  P2(A);

  mtl::dense_vector<double> x(N, 1.0), b(N);

  b = A * x;

  for (int l= 2; l <= 10; l++) {
      x= 0;
      itl::cyclic_iteration<double> iter(b, 100, 1.e-6, 0.0, 100);
      boost::timer t;
      bicgstab_ell(A, x, b, P, P2, iter, l);
      std::cout << "BiCGStab(" << l << ") took " << t.elapsed() << "s.\n";
  }
 
  return 0;
}
