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

	// Inner domain
	for (int i= 1; i < m-1; i++)
	    for (int j= 1, k= i * n + j; j < n-1; j++, k++) 
		w[k]= 4 * v[k] - v[k-n] - v[k+n] - v[k-1] - v[k+1]; 
	    
	// Upper border
	for (int j= 1; j < n-1; j++) 
	    w[j]= 4 * v[j] - v[j+n] - v[j-1] - v[j+1];

	// Lower border
	for (int j= 1, k= (m-1) * n + j; j < n-1; j++, k++) 
	    w[k]= 4 * v[k] - v[k-n] - v[k-1] - v[k+1]; 
	
	// Left border
	for (int i= 1, k= n; i < m-1; i++, k+= n)
	    w[k]= 4 * v[k] - v[k-n] - v[k+n] - v[k+1]; 

	// Right border
	for (int i= 1, k= n+n-1; i < m-1; i++, k+= n)
	    w[k]= 4 * v[k] - v[k-n] - v[k+n] - v[k-1]; 

	// Corners
	w[0]= 4 * v[0] - v[1] - v[n];
	w[n-1]= 4 * v[n-1] - v[n-2] - v[2*n - 1];
	w[(m-1)*n]= 4 * v[(m-1)*n] - v[(m-2)*n] - v[(m-1)*n+1];
	w[m*n-1]= 4 * v[m*n-1] - v[m*n-2] - v[m*n-n-1];

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
    typedef mtl::dense_vector<double> vt;
    
    mtl::compressed2D<double> A0;
    laplacian_setup(A0, 4, 5);
    cout << "A0 is\n" << A0 << endl;

    vt v(20);
    iota(v);
    cout << "v is " << v << endl;

    vt w1(A0 * v);
    cout << "A0 * v is " << w1 << endl;

    poisson2D_dirichlet A(4, 5);
    vt  w2;
    w2= A * v;
    cout << "A * v is " << w2 << endl;

    if (one_norm(vt(w1 - w2)) > 0.001) throw "Wrong result";


    return 0;
}
