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
 
using namespace std;  

template <typename Matrix, typename Vector>
void inline power_iteration(const Matrix& A, Vector& v, double tau)
{
    assert(num_rows(A) == num_cols(A));               // A should be square    
    v*= 1. / two_norm(v);                             // Normalize v
    Vector v2(size(v));
    do {	 
	swap(v, v2);                                  // Keep old value in v2    
	v= A * v2;           	
	v*= 1. / two_norm(v);                         // Normalize 
    } while (two_norm(Vector(v - v2)) >= tau);
}

template <typename Matrix, typename Vector>
bool inline check_eigenvector(const Matrix& A, Vector& v, double tau)
{    
    Vector w(A * v);
    std::cout << "A * v is " << w << endl;

    typename mtl::Collection<Vector>::value_type alpha= two_norm(w) / two_norm(v);
    std::cout << "Eigenvalue alpha for v is " << alpha << endl;
    Vector w2(alpha * v);
    std::cout << "alpha * v is " << w2 << endl;

    bool close= two_norm(Vector(w - w2)) / two_norm(Vector(w + w2)) < tau;
    std::cout << "The results are " << (close ? "similar.\n" : "different.\n") << endl;
    return close;
}

int main(int, char**)
{
    double a_value[4][4] = {{0, 0, 1, .5},
			    {1/3., 0, 0, 0},
			    {1/3., .5, 0, .5},
			    {1/3., .5, 0, 0}};
    mtl::dense2D<double> A(a_value);
    mtl::dense_vector<double> v(4, 1.0);

    power_iteration(A, v, 0.00001);
    check_eigenvector(A, v, 0.001);

    mtl::compressed2D<double> B(9, 9);
    {
	mtl::mat::inserter<mtl::compressed2D<double> > ins(B);
	ins[0][1] << .2; ins[0][4] << .5;
	ins[1][0] << .5; ins[1][3] << 1; ins[1][4] << .5; ins[1][5] << .25; ins[1][6] << 1/3.;
	ins[2][2] << .1; ins[2][5] << .25;
	ins[3][1] << .2;
	ins[4][0] << .5; ins[4][1] << .2;
	ins[5][1] << .2; ins[5][6] << 1/3.; ins[5][8] << 1; 
	ins[6][1] << .2; ins[6][5] << .25; ins[6][7] << 1; 
	ins[7][6] << 1/3.; 
	ins[8][5] << .25; 
    }
    mtl::dense_vector<double> w(9, 1.0);

    power_iteration(B, w, 0.00001);
    check_eigenvector(B, w, 0.001);

    return 0;
}
