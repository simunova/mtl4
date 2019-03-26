// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschrÃ¤nkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/timer.hpp>


using namespace std;

template <typename Matrix, typename Value, typename Vector>
void test_vector(const Matrix& A, const Value& alpha, const Vector& v, double tol, int i)
{
    Vector diff(A*v-alpha*v), v1(A*v), v2(alpha*v); // , diff(v1-v2);
    if (size(v1) < 17) 	
	cout << "A*v is     " << v1 << "\nalpha*v is " << v2 << '\n';
    if (two_norm(diff) > tol) cout << "two_norm(difference) of the " << i << "-th eigenvector is " << two_norm(diff) << ", two_norm(A*v) is " << two_norm(v1) << ", two_norm(alpha*v) is " << two_norm(v2) << '\n'; // throw "wrong eigenvector";
}


int main(int argc, char** argv) 
{
    using namespace mtl;

    int select= 1, sub= 600;
    if (argc > 1)
	select= atoi(argv[1]);
    assert(select >= 1 && select <= 2);

    if (argc > 2)
	sub= atoi(argv[2]);
    
    string fname= string("../../../../../branches/data/matrix_market/Partha") + char('0' + select) + ".mtx";
    

    dense2D<double>    A0(io::matrix_market(fname.c_str())), A(clone(A0[irange(sub)][irange(sub)]));
    //    cout << "Size of A is " << num_rows(A) << " x " << num_cols(A) << '\n';
   
    boost::timer tri_time;
    dense2D<double>    C(hessenberg_factors(A)), D(clone(bands(C, -1, 2))), Q(num_rows(D), num_rows(D));
    cout << "The tridiagonal matrix is\n" << D[irange(10)][irange(10)] << "This took " << tri_time.elapsed() << "s.\n";
   
    tri_time.restart();
    dense_vector<double>       lambda(num_rows(D));

    cuppen(D, Q, lambda);
    cout << "Q is\n" << Q[irange(10)][irange(10)] << "This took " << tri_time.elapsed() << "s.\n";
    // std::cout << "The eigenvalues are " << lambda << "\n";
#if 0  
    for (unsigned i= 0; i < num_rows(D); i++)
	test_vector(D, lambda[i], mtl::dense_vector<double>(Q[mtl::iall][i]), 1e-4, i);
#endif
    return 0;
}
