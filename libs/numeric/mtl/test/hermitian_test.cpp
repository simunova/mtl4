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
#include <typeinfo>
#include <boost/numeric/mtl/mtl.hpp>

//#include <boost/numeric/mtl/matrix/hermitian_view.hpp>

using namespace std;

template <typename Matrix>
void test_type(Matrix)
{
    
    cout << "Type " << typeid(Matrix).name() 
	 << "\n     has ashape " << typeid(typename mtl::ashape::ashape<Matrix>::type).name() << '\n';
    cout << "Type " << typeid(const Matrix).name() 
	 << "\n     has ashape " << typeid(typename mtl::ashape::ashape<const Matrix>::type).name() << '\n';
}

template <typename Matrix>
void test(const Matrix&, const char* name)
{
    cout << name 
	    // << " ... " << typeid(typename mtl::OrientedCollection<mtl::mat::hermitian_view<Matrix> >::orientation).name() 
	 << '\n';

    cout << name << '\n';
    Matrix A(3, 3), B(3, 3), C(3, 3);
    // Stupid test with real values --> use complex
    A= 2, 3, 4,
       1, 2, 3,
       3, 4, 5;
    B= 3, 4, 5,
       2, 3, 4,
       1, 2, 3;

    C= hermitian(B);
    cout << "hermitian(B) is\n" << C;

    C= A * hermitian(B);
    cout << "A * hermitian(B) is\n" << C;

    typedef typename mtl::Collection<Matrix>::value_type value_type;
    mtl::dense_vector<value_type> v(3), w(3, 1.0);

#if 0 // understanding ashape of const views
    const Matrix K(A);
    test_type(K);
    test_type(trans(K));
    test_type(conj(K));
    test_type(conj(trans(K)));
    test_type(trans(conj(K)));
#endif

    v= hermitian(A) * w;
    v= mtl::mat::conj(trans(A)) * w;
    v= trans(mtl::mat::conj(A)) * w;
    v= trans(A) * w;
}


int main(int, char**)
{
    using namespace mtl;
    dense2D<double>                                      dr;
    dense2D<double, mat::parameters<col_major> >      dc;
    morton_dense<double, recursion::morton_z_mask>       mzd;
    morton_dense<double, recursion::doppled_2_row_mask>  d2r;
    compressed2D<double>                                 cr;
    compressed2D<double, mat::parameters<col_major> > cc;

    dense2D<complex<double> >                            drc;
    compressed2D<complex<double> >                       crc;

    test(dr, "Dense row major");
#if 0
    test(dc, "Dense column major");
    test(mzd, "Morton Z-order");
    test(d2r, "Hybrid 2 row-major");
    test(cr, "Compressed row major");
    test(drc, "Dense row major complex");

    // For better readability I don't want finish with a complex
    test(cc, "Compressed column major");
#endif
    return 0;
}
