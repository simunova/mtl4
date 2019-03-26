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

void test(int n)
{
    mtl::compressed2D<double> A(n*n, n*n);
    {
        mtl::mat::inserter<mtl::compressed2D<double> > ins(A, 5);
	for (int i= 0; i < n; i++) {
	    int r= i * n, re= r+n;
	    if (i > 0)
		for (int j= r; j < re; j++)
		    ins[j][j-n] << -1.0;
	    ins[r][r] << 4.0;
	    ins[r][r+1] << -1;
	    for (int j= r+1; j < re-1; j++) {
		ins[j][j-1] << -1.0;
		ins[j][j] << 4.0;
		ins[j][j+1] << -1.0;
	    }
	    ins[re-1][re-2] << -1.0;
	    ins[re-1][re-1] << 4.0;
	    if (i < n-1)
		for (int j= r; j < re; j++)
		    ins[j][j+n] << -1.0;
	}
    }
    if (n < 9) cout << "A is \n" << A;
    MTL_THROW_IF(A.nnz() != unsigned(5*n*n - 4*n), mtl::runtime_error("Wrong number of non-zeros."));
    for (int i= 0, ns= n*n; i < ns; i++) {
	MTL_THROW_IF(A[i][i] != 4.0, mtl::runtime_error("Wrong value"));
	MTL_THROW_IF(i%n != 0 && A[i][i-1] != -1.0, mtl::runtime_error("Wrong value"));
	MTL_THROW_IF((i+1)%n != 0 && A[i][i+1] != -1.0, mtl::runtime_error("Wrong value"));
	MTL_THROW_IF(i-n >= 0 && A[i][i-n] != -1.0, mtl::runtime_error("Wrong value"));
	MTL_THROW_IF(i+n < ns && A[i][i+n] != -1.0, mtl::runtime_error("Wrong value"));
    }
	
}

 
int main(int, char**)
{
    test(9);

    return 0;
}
