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
#include <boost/timer.hpp>
#include <boost/numeric/mtl/mtl.hpp>


using namespace std;

template <typename Matrix, typename Vector>
typename mtl::Collection<Matrix>::value_type
matgen (Matrix& A, Vector& b)
{
    size_t                                        n= num_rows(A);
    typename mtl::Collection<Matrix>::value_type  norma= 0, tmp;
    int                                           init= 1325;

    for (size_t i = 0; i < n; i++) 
	for (size_t j = 0; j < n; j++) {
	    init = 3125*init % 65536;
	    A[j][i]= tmp= (init - 32768.0) / 16384.0;
	    norma = tmp > norma ? tmp : norma;
	}
    for (size_t i = 0; i < n; i++) 
	b[i]= mtl::sum(A[i][mtl::iall]);

    return norma;
}




template <typename Matrix>
void test(Matrix& A, const char* name)
{
    cout << "\n" << name << "\n";

    typedef typename mtl::Collection<Matrix>::value_type  Scalar;
    typedef typename mtl::dense_vector<Scalar>            Vector;

    Vector b(num_rows(A), 0.0), x;
    matgen(A, b);

    boost::timer s;
    x= lu_solve(A, b);
    double t= s.elapsed(), n= num_rows(A), ops= 2./3.*n*n*n + 2*n*n;
    std::cout << "LU solution took " << t << "s. This corresponds to "
	      << ops / t / 1e6 << " MFlops.\n";

#if 0
    cout << "A is\n" << A << '\n';
    cout << "b is " << b << '\n';
    cout << "x is " << x << '\n';

    b= A * x;
    cout << "A*x is " << b << '\n';
#endif
}



int main(int, char**)
{
    using namespace mtl;
    unsigned size= 1000;
    
    dense2D<double>                                      dr(size, size);
    dense2D<complex<double> >                            dz(size, size);
    dense2D<double, mat::parameters<col_major> >      dc(size, size);

    test(dr, "Row-major dense");
    //test(dz, "Row-major dense with complex numbers");
    test(dc, "Column-major dense");

    return 0;
}
