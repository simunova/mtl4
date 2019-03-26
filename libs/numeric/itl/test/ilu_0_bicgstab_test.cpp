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

template <typename At, typename Lt, typename Ut>
void dense_ilu_0(const At& As, const Lt& Ls, const Ut& Us)
{
    mtl::dense2D<double> LU(As);
     
    const std::size_t n= num_rows(LU);
    for (std::size_t i= 1; i < n; i++) 
	for (std::size_t k= 0; k < i; k++) {
	    LU[i][k]/= LU[k][k];
	    for (std::size_t j= k + 1; j < n; j++)
		if (LU[i][j] != 0)
		    LU[i][j]-= LU[i][k] * LU[k][j];
	}
    std::cout << "Factorizing A = \n" << As << "-> LU = \n" << LU;
    // std::cout << "L = \n" << Ls << "\nU = \n" << Us;

    if (std::abs(LU[3][2] - Ls[3][2]) > 0.001) throw "Wrong value in L for sparse ILU(0) factorization";

    if (std::abs(LU[3][3] - 1. / Us[3][3]) > 0.001) throw "Wrong value in U for sparse ILU(0) factorization";
}


int main()
{
    // For a more realistic example set sz to 1000 or larger
    const int size = 3, N = size * size; 

    typedef mtl::compressed2D<double>  matrix_type;
    mtl::compressed2D<double>          A(N, N), dia(N, N);
    laplacian_setup(A, size, size);
    // dia= 1.0; A+= dia;
    
   
    itl::pc::ilu_0<matrix_type>        P(A);
    mtl::dense_vector<double>          x(N, 1.0), b(N);
    
    if(size > 1 && size < 4)
	dense_ilu_0(A, P.get_L(), P.get_U());

    b = A * x;
    x= 0;
    
    itl::cyclic_iteration<double> iter(b, N, 1.e-6, 0.0, 5);
    bicgstab(A, x, b, P, iter);
    
    return 0;
}
