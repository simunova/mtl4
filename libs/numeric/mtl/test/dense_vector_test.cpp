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
#include <cmath>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/operation/operators.hpp>
#include <boost/numeric/mtl/operation/dot.hpp>
#include <boost/numeric/mtl/operation/pow.hpp>


using namespace std;  


template <typename Vector>
void one_d_iteration(char const* name, Vector & vector, size_t check_index, typename Vector::value_type check)
{
    namespace traits = mtl::traits;
    typename traits::index<Vector>::type                               index(vector);
    typename traits::value<Vector>::type                               value(vector); 
    typedef  mtl::tag::nz                                              tag;
    typedef typename traits::range_generator<tag, Vector>::type        cursor_type;
    typedef typename traits::range_generator<tag, Vector>::complexity  complexity;

    cout << name << "\nElements: " << complexity() << '\n';
    for (cursor_type cursor = mtl::begin<tag>(vector), cend = mtl::end<tag>(vector); cursor != cend; ++cursor) {
	// cout << "vector[" << index(*cursor) << "] = " << value(*cursor) << '\n';
	if (index(*cursor) == check_index && value(*cursor) != check) 
	    throw "wrong check value";
    }
}
    

template <typename VectorU, typename VectorV, typename VectorW>
void test(VectorU& u, VectorV& v, VectorW& w, const char* name)
{
    u= 3.0; v= 4.0; w= 5.0;

    std::cout << "\n\n";
    one_d_iteration(name, u, 2, (typename VectorU::value_type)(3.0));

    std::cout << name << "  --- u= v + w:" << std::endl;
    u= v + w;
    cout << "u: " << u << std::endl;
    MTL_THROW_IF(u[0] != 9.0, mtl::runtime_error("wrong"));

    std::cout << name << "  --- u= ele_prod(v, w):" << std::endl;
    u= ele_prod(v, w);
    cout << "u: " << u << std::endl;
    MTL_THROW_IF(u[0] != 20.0, mtl::runtime_error("wrong"));

    std::cout << name << "  --- u= ele_prod(v, w) + w:" << std::endl;
    u= ele_prod(v, w) + w;
    cout << "u: " << u << std::endl;
    MTL_THROW_IF(u[0] != 25.0, mtl::runtime_error("wrong"));

    std::cout << name << "  --- u= ele_quot(v, w):" << std::endl;
    u= ele_quot(v, w);
    cout << "u: " << u << std::endl;
    MTL_THROW_IF(abs(u[0] - 0.8) > 0.01, mtl::runtime_error("wrong"));

    std::cout << name << "  --- u= ele_quot(v+w, w):" << std::endl;
    u= ele_quot(v+w, w);
    cout << "u: " << u << std::endl;
    MTL_THROW_IF(abs(u[0] - 1.8) > 0.01, mtl::runtime_error("wrong"));


    std::cout << name << "  --- u= v + w + v + w:" << std::endl;
    u= v + w + v + w;
    cout << "u: " << u << std::endl;
    MTL_THROW_IF(u[0] != 18.0, mtl::runtime_error("wrong"));

    std::cout << name << "  --- u= w + (v= w + w);:" << std::endl;
    u= w + (v= w + w);
    cout << "u: " << u << "v: " << v << std::endl;
    MTL_THROW_IF(v[0] != 10.0, mtl::runtime_error("v wrong"));
    MTL_THROW_IF(u[0] != 15.0, mtl::runtime_error("u wrong"));

    std::cout << name << "  --- u= (v= w + w) + v;:" << std::endl;
    v= 4.0; w= 5.0; u= (v= w + w) + v;
    cout << "u: " << u << "v: " << v << std::endl;
    MTL_THROW_IF(v[0] != 10.0, mtl::runtime_error("v wrong"));
    MTL_THROW_IF(u[0] != 20.0, mtl::runtime_error("u wrong"));

    std::cout << name << "  --- w= 4; u-= (v= w + w) - w;:" << std::endl;
    w= 4; u-= (v= w + w) - w;
    cout << "u: " << u << "v: " << v << std::endl;
    MTL_THROW_IF(v[0] != 8.0, mtl::runtime_error("v wrong"));
    MTL_THROW_IF(u[0] != 16.0, mtl::runtime_error("u wrong")); // for -=

    
    std::cout << name << "  --- v= 3*u + 4*w;:" << std::endl;
    v= 3*u + 4*w;
    cout << "v: " << v << std::endl;
    MTL_THROW_IF(v[0] != 64.0, mtl::runtime_error("v wrong"));

    std::cout << name << "  --- u= 3; v= 4; u+= v+= 3.0 * (w= 5);:" << std::endl;
    u= 3; v= 4; 
    u+= v+= 3.0 * (w= 5.0);
    cout << "u: " << u << "v: " << v << std::endl;
    MTL_THROW_IF(v[0] != 19.0, mtl::runtime_error("v wrong"));
    MTL_THROW_IF(u[0] != 22.0, mtl::runtime_error("u wrong"));

    std::cout << name << "  --- u= 3; v= 4; w=5; u+= (v*= 3.0) + (w*= 2.0);:" << std::endl;
    u= 3; v= 4; w=5; u+= (v*= 3.0) + (w*= 2.0);
    cout << "u: " << u << "v: " << v << "w: " << w << std::endl;
    MTL_THROW_IF(v[0] != 12.0, mtl::runtime_error("v wrong"));
    MTL_THROW_IF(w[0] != 10.0, mtl::runtime_error("v wrong"));
    MTL_THROW_IF(u[0] != 25.0, mtl::runtime_error("u wrong"));

    std::cout << name << "  --- u= 3; v= 4; w=5; u+= dot(v, w) * w;:" << std::endl;
    u= 3; v= 4; w=5; u+= dot(v, w) * w;
    cout << "u: " << u << "v: " << v << "w: " << w << std::endl;
    MTL_THROW_IF(u[0] != 503.0, mtl::runtime_error("u wrong"));

    std::cout << name << "  --- u+= dot<12>(v, w) * w;:" << std::endl;
    u+= mtl::dot<12>(v, w) * w;
    cout << "u: " << u << "v: " << v << "w: " << w << std::endl;
    MTL_THROW_IF(u[0] != 1003.0, mtl::runtime_error("u wrong"));
    
    std::cout << name << "  --- u= pow(v, 2.0); // with v= 1, 2, 3, 4, 5;" << std::endl;
    v= 1, 2, 3, 4, 5;
    u= pow(v, 2.0);
    cout << "u: " << u << "v: " << v << std::endl;
    MTL_THROW_IF(abs(u[4] - 25) > 0.0001, mtl::runtime_error("u wrong")); 
}
 

int main(int , char**)
{
    using mtl::vec::parameters;
    using namespace mtl;

    dense_vector<float>   u(5), v(5), w(5);
    dense_vector<double>  x(5), y(5), z(5);
    dense_vector<std::complex<double> >  xc(5), yc(5), zc(5);

    std::cout << "Testing vector operations\n";

    test(u, v, w, "test float");
    test(x, y, z, "test double");
    test(u, x, y, "test float, double mixed");
#if 0
    // test is not designed for complex
    test(xc, yc, zc, "test complex<double>");
    test(x, yc, zc, "test complex<double>, double mixed");
#endif

    dense_vector<float, parameters<row_major> >   ur(5), vr(5), wr(5);
    test(ur, vr, wr, "test float in row vector");
    
    // test(ur, v, wr, "test float in mixed vector (shouldn't work)"); 

    return 0;
}
 














