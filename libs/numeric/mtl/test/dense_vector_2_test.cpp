/*
 *  vector_right_scale_test.cpp
 *  MTL
 *
 *  Created by Hui Li (huil@Princeton.EDU)
 *
 */

#include <iostream>
#include <cmath>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/operation/operators.hpp>
#include <boost/numeric/mtl/operation/dot.hpp>


using namespace std;  


template <typename VectorU, typename VectorV, typename VectorW>
void test(VectorU& u, VectorV& v, VectorW& w, const char* name)
{
    using mtl::dot;

	// test right scale
    u= 3.0; v= 4.0; w= 5.0;
	
    std::cout << name << "  --- v= u*3 + w*4;:\n"; std::cout.flush();
    v= u*3 + w*4;
    cout << "v: " << v << "\n"; std::cout.flush();
    MTL_THROW_IF(v[0] != 29.0, mtl::runtime_error("v wrong"));
	
    std::cout << name << "  --- u= 3; v= 4; u+= v+= (w= 5) * 3.0;:\n"; std::cout.flush();
    u= 3; v= 4; 
    u+= v+= (w= 5.0) * 3.0;
    cout << "u: " << u << "v: " << v << "\n"; std::cout.flush();
    MTL_THROW_IF(v[0] != 19.0, mtl::runtime_error("v wrong"));
    MTL_THROW_IF(u[0] != 22.0, mtl::runtime_error("u wrong"));
	
    std::cout << name << "  --- u= 3; v= 4; w=5; u+= w * dot(v, w);:\n"; std::cout.flush();
    u= 3; v= 4; w=5; u+= w * dot(v, w);
    cout << "u: " << u << "v: " << v << "w: " << w << "\n"; std::cout.flush();
    MTL_THROW_IF(u[0] != 503.0, mtl::runtime_error("u wrong"));
	
    std::cout << name << "  --- u+= w * dot<12>(v, w);:\n"; std::cout.flush();
    u+= w * dot<12>(v, w);
    cout << "u: " << u << "v: " << v << "w: " << w << "\n"; std::cout.flush();
    MTL_THROW_IF(u[0] != 1003.0, mtl::runtime_error("u wrong"));

	// test divide by scalar
    u= 3.0; v= 4.0; w= 5.0;
	
    std::cout << name << "  --- v= u/3 + w/5;:\n"; std::cout.flush();
    v= u/3 + w/5;
    cout << "v: " << v << "\n"; std::cout.flush();
    MTL_THROW_IF(v[0] != 2.0, mtl::runtime_error("v wrong"));
	
    std::cout << name << "  --- u= 3; v= 4; u+= v+= (w= 6.0) / 3.0;:\n"; std::cout.flush();
    u= 3; v= 4; 
    u+= v+= (w= 6.0) / 3.0;
    cout << "u: " << u << "v: " << v << "\n"; std::cout.flush();
    MTL_THROW_IF(v[0] != 6.0, mtl::runtime_error("v wrong"));
    MTL_THROW_IF(u[0] != 9.0, mtl::runtime_error("u wrong"));
	
}

int main(int, char**)
{
    using mtl::vec::parameters; using mtl::dense_vector;
	
    dense_vector<float>   u(5), v(5), w(5);
    dense_vector<double>  x(5), y(5), z(5);
    dense_vector<std::complex<double> >  xc(5), yc(5), zc(5);
	
    std::cout << "Testing vector right scaling operations\n";
	
    test(u, v, w, "test float");
    test(x, y, z, "test double");
    test(u, x, y, "test float, double mixed");
#if 0
    // test is not designed for complex
    test(xc, yc, zc, "test complex<double>");
    test(x, yc, zc, "test complex<double>, double mixed");
#endif
	
    dense_vector<float, parameters<mtl::row_major> >   ur(5), vr(5), wr(5);
    test(ur, vr, wr, "test float in row vector");
    
    // test(ur, v, wr, "test float in mixed vector (shouldn't work)"); 
	
    return 0;
}
