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
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

struct poisson2D_dirichlet
{
    poisson2D_dirichlet(int m, int n) : m(m), n(n) {}

    template <typename Vector>
    Vector operator*(const Vector& v) const
    {
	assert(int(size(v)) == m * n);
	Vector w(m * n);
	
	for (int i= 0; i < m; i++)
	    for (int j= 0; j < n; j++) {
		int k= i * n + j; // offset
		w[k]= 4 * v[k];
		if (i > 0) w[k]-= v[k-n];   // upper neighbor
		if (i < m-1) w[k]-= v[k+n]; // lower neighbor
		if (j > 0) w[k]-= v[k-1];   // left neighbor
		if (j < n-1) w[k]-= v[k+1]; // right neighbor
	    }
	return w;
    }
    int m, n;
};

namespace mtl { namespace ashape {
    template <> struct ashape_aux<poisson2D_dirichlet> 
    {	typedef nonscal type;    };
}}


int main()
{
    mtl::vampir_trace<9999> tracer;
  // For a more realistic example set size to 1000 or larger
  const int size = 10, N = size * size;

  poisson2D_dirichlet A(size, size);

  itl::pc::identity<poisson2D_dirichlet>     P(A);

  mtl::dense_vector<double> x(N, 1.0), b(N);

  b = A * x;
  x= 0;
  itl::cyclic_iteration<double> iter(b, 100, 1.e-11, 0.0, 5);
  cg(A, x, b, P, iter);

  return 0;
}
