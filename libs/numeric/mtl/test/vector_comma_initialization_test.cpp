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

template <typename Vector>
void test(Vector& v, const char* name)
{
    std::cout << name << "\n";

    using mtl::Collection;    
    typedef typename Collection<Vector>::value_type   value_type;
    
    v= 7.;
    std::cout << "After v= 7; v == " << v << "\n";
    MTL_THROW_IF(v[2] != value_type(7.), mtl::runtime_error("Constant assignment wrong"));

    v= 8., 45u, -3;
    std::cout << "After v= 8., 45., -3.; v == " << v << "\n";
    MTL_THROW_IF(v[2] != value_type(-3), mtl::runtime_error("Comma assignment wrong"));

    v= 6.;
    std::cout << "After v= 6; v == " << v << "\n";
    MTL_THROW_IF(v[2] != value_type(6.), mtl::runtime_error("Constant reassignment wrong"));

    Vector w(3);
    v= w= 5.;
    std::cout << "After v= w= 5; v == " << v << ", w == " << w << "\n\n";
    MTL_THROW_IF(w[2] != value_type(5.), mtl::runtime_error("Constant assignment (of w) wrong"));
    MTL_THROW_IF(v[2] != value_type(5.), mtl::runtime_error("Constant reassignment wrong"));

#if 0
    v+= 8., 45u, -3;
    std::cout << "After v+= 8., 45., -3.; v == " << v << "\n";
    if (v[0] != value_type(13))
	throw "Comma assignment wrong";
    if (v[2] != value_type(2))
	throw "Comma assignment wrong";
#endif

}



int main(int, char**)
{
    mtl::dense_vector<float>                  v(3);
    mtl::dense_vector<std::complex<float> >   w(3);

    test(v, "dense_vector<float>");
    test(w, "dense_vector<std::complex<float> >");

    return 0;
}
