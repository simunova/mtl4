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

typedef std::complex<double> cmp;

template <typename Vector>
void test(const Vector&, const char* name)
{
    mtl::multi_vector<Vector> A(5, 5);
    A= 3.0;
    A[3][2]= cmp(0, 1);
    cout << name << ":\n A after initialization is \n" << A << "\n";

    mtl::dense_vector<cmp>  b(5, 4.0), x;
    x= trans(A) * b;

    cout << "x = " << x << "\n";
    MTL_THROW_IF(x[1] != 12.0, mtl::runtime_error("Wrong value"));
    MTL_THROW_IF(x[2] != cmp(12.0, 4.0), mtl::runtime_error("Wrong value"));

    mtl::dense_vector<cmp, mtl::vec::parameters<mtl::row_major> >  br(5, 4.0), xr;
    xr= br * A;

    cout << "xr = " << xr << "\n";
    MTL_THROW_IF(xr[1] != 12.0, mtl::runtime_error("Wrong value"));
    MTL_THROW_IF(xr[2] != cmp(12.0, 4.0), mtl::runtime_error("Wrong value"));

    xr= br * hermitian(A);  // corresponds to trans(x)= trans( conj(A) * b )
    cout << "xr = " << xr << "\n";
    MTL_THROW_IF(xr[2] != 12.0, mtl::runtime_error("Wrong value"));
    MTL_THROW_IF(xr[3] != cmp(12.0, -4.0), mtl::runtime_error("Wrong value"));

}


int main(int, char**)
{
    using namespace mtl;

    dense_vector<cmp>    v;
    test(v, "dense_vector<double>");

    return 0;
}
