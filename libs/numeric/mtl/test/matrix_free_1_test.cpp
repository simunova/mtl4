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

#include <typeinfo>
#include <iostream>
#include <cassert>
#include <boost/numeric/mtl/mtl.hpp>

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

int main(int, char**)
{
    using namespace std;
    
    mtl::compressed2D<double> A0;
    laplacian_setup(A0, 4, 5);
    cout << "A0 is\n" << A0 << endl;

    mtl::dense_vector<double> v(20);
    iota(v);
    cout << "v is " << v << endl;

    mtl::dense_vector<double> w1(A0 * v);
    cout << "A0 * v is " << w1 << endl;

    poisson2D_dirichlet A(4, 5);
    mtl::dense_vector<double> w2(A * v);
    cout << "A * v is " << w2 << endl;

    w1-= w2;
    if (one_norm(w1) > 0.001) throw "Wrong result";

    return 0;
}
