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

typedef std::complex<double>      cdouble;

template <typename Matrix>
void test(Matrix& A, const char* name)
{
    //using mtl::conj;
    const unsigned                    xd= 2, yd= 5, n= xd * yd;
    A.change_dim(n, n);
    laplacian_setup(A, xd, yd); 

    A*= cdouble(1, -1);
    std::cout << name << "\nconj(A) is\n" << with_format(mtl::conj(A), 7, 1) << "\n";

    mtl::dense_vector<cdouble> x(n),Ax(n);
    x=cdouble(1,2);
    
    // Ax= mtl::mat::conj(A) * x;
    Ax= mtl::conj(A) * x;
    std::cout << "conj(A) * x is " << Ax << "\n";
    
    Ax=trans(A) * x;
    std::cout << "trans(A) * x is " << Ax << "\n";

    Ax=hermitian(A) * x;
    std::cout << "hermitian(A) * x is " << Ax << "\n";
}


int main(int, char**)
{
    mtl::compressed2D<cdouble>             crc;

    test(crc, "Compressed row major complex");

    return 0;
}
