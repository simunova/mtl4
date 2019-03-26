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
#include <complex>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;  

#if defined(MTL_WITH_INITLIST) && defined(MTL_WITH_AUTO) && defined(MTL_WITH_RANGEDFOR)
#include <initializer_list>

template <typename Matrix>
void test(const char* name)
{
    const Matrix A={{3, 4}, 
		    {5, 6}};

    cout << "Test " << name << ", A(initialization):\n" << A << "\n"; 
    MTL_THROW_IF(A[1][0] != 5.0, mtl::runtime_error("wrong value or position"));

    Matrix B;
    B= {{3, 4}, 
	{5, 6}};
    cout << "B(assignment):\n" << B << "\n"; 
    MTL_THROW_IF(B[1][0] != 5.0, mtl::runtime_error("wrong value or position"));
}
#else
template <typename Matrix>
void test(const char* ) {}
#endif 

int main(int , char**)
{
    using mtl::mat::parameters;
    using namespace mtl;
    
    test<dense2D<double> >("dense2D<double>");
    test<dense2D<float> >("dense2D<float>");
    test<dense2D<complex<double> > >("dense2D<complex<double>>");
    test<dense2D<float, mat::parameters<col_major> > >("dense2D<float> col_major");

    test<compressed2D<float> >("compressed2D<float>");
    test<compressed2D<float, mat::parameters<col_major> > >("compressed2D<float> col_major");
    
    test<morton_dense<double, recursion::morton_z_mask> >("Morton Z-order");

    return 0;
}
 














