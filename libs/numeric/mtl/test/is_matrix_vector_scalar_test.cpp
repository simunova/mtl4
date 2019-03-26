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
#include <complex>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;
    
template <typename T>
void test(int expected)
{
    cout << typeid(T).name() << " is " << (mtl::traits::is_scalar<T>::value ? "" : "not ") << "scalar\n";
    cout << "Is " << (mtl::traits::is_vector<T>::value ? "" : "not ") << "vector\n";
    cout << "Is " << (mtl::traits::is_matrix<T>::value ? "" : "not ") << "matrix\n\n";

    MTL_THROW_IF((expected != int(mtl::traits::is_scalar<T>::value) + 2 * int(mtl::traits::is_vector<T>::value)
		  + 3 * int(mtl::traits::is_matrix<T>::value)), mtl::unexpected_result());
}

int main(int, char**)
{
    using namespace mtl;
    
    test<int>(1);
    test<float>(1);
    test<double>(1);
    test<std::complex<double> >(1);
    
    test<dense_vector<double> >(2);
    
    test<dense2D<double> >(3);
    test<compressed2D<double> >(3);

    return 0;
}
