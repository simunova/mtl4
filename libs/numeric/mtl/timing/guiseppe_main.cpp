#include <iostream>
#include <stdio.h>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/timer.hpp>

#include "ScalarVariations.hpp"

/*
  First run: 170us dynamic types (220us old value)
              52us static types (r6809)
	      56us static types with unrolled by hand (r6810)
	      93us without my unrolled matrix product
	      55us with my unrolled matrix vector product (r6813)
*/

using namespace mtl;
using namespace mtl::matrix;
using namespace mtl::vector;

#define STATIC_TYPES

#ifdef STATIC_TYPES
   typedef dense_vector<double, parameters<tag::col_major, fixed::dimension<3>, true> > vec;
   typedef mat::parameters<tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<3, 3>, true> mat_para;
   #define VEC_ARG
#else
   typedef dense_vector<double> vec;
   typedef mat::parameters<> mat_para;
   #define VEC_ARG (3)
#endif
   typedef dense2D<double, mat_para> mat;
using namespace std;

int main()
{


    vec X VEC_ARG; X[0]=21; X[1]=2;  X[2]=31;

    vec F1 VEC_ARG; F1[0]=21.; F1[1]=2.;  F1[2]=31.;
    vec F2 VEC_ARG; F2[0]=31.; F2[1]=0.;  F2[2]=32.;
    vec F3 VEC_ARG; F3[0]=41.; F3[1]=42.; F3[2]=33.;
    vec F4 VEC_ARG; F4[0]=21.; F4[1]=24.; F4[2]=43.;

    vec A1 VEC_ARG; A1[0]=11.; A1[1]=25.; A1[2]=6.;
    vec A2 VEC_ARG; A2[0]=12.; A2[1]=26.; A2[2]=7.;
    vec A3 VEC_ARG; A3[0]=13.; A3[1]=27.; A3[2]=8.;
    vec A4 VEC_ARG; A4[0]=14.; A4[1]=28.; A4[2]=9.;

    vec U1 VEC_ARG; U1[0]=11.; U1[1]=25.; U1[2]=6.;
    vec U2 VEC_ARG; U2[0]=12.; U2[1]=26.; U2[2]=7.;
    vec U3 VEC_ARG; U3[0]=13.; U3[1]=27.; U3[2]=8.;
    vec U4 VEC_ARG; U4[0]=14.; U4[1]=28.; U4[2]=9.;

    Variations<mat,vec> vars;  // solo qua

    mat q1 = vars.Q1(A1);
    std::cout << q1 << std::endl;

    mat q2 = vars.Q1(A2);
    std::cout << (q2) << std::endl;

    mat q12 = vars.Q2(A1,A2);
    std::cout << q12 << std::endl;

    mat q123 = vars.Q3(A1,A2,A3);
    std::cout << q123 << std::endl;

    mat q1234 = vars.Q4(A1,A2,A3,A4);
    std::cout << q1234 << std::endl;

    mat c = vars.skew(q1234);
    std::cout << c << std::endl;

    mat w_1 = vars.w1(A1,F1);
    std::cout << w_1 << std::endl;

    mat w_12 = vars.w2(A1,F1,A2,F2);
    std::cout <<  w_12 << std::endl;

    mat w_123 = vars.w3(A1,F1,A1,F1,A1,F1);
    std::cout <<  w_123 << std::endl;

    mat w_1234 = vars.w4(A1,F1,A1,F1,A3,F3,A4,F4);
    std::cout << w_1234 << std::endl;

    vec _f1 = vars.f1(A1,F1);
    std::cout << _f1 << std::endl;

    vec _f12 = vars.f2(A1,F1,A2,F2);
    std::cout << _f12 << std::endl;

    std::cout << vars.c2(A1,F1,A1,F1) << std::endl;
    std::cout << vars.c2(A1,F1,A2,F2) << std::endl;
    std::cout << vars.c2(A2,F2,A3,F3) << std::endl;

    vec _f123 = vars.f3(A1,F1,A2,F2,A1,F1);
    std::cout << _f123 << std::endl;

    vec _f1234 = vars.f4(A1,F1,A2,F2,A1,F1,A4,F4);
    std::cout << "f4=" << _f1234 << std::endl;

    vec _u1 = vars.u1(X, A1, U1);
    std::cout << "u1=" << _u1 << std::endl;

    vec _u12 = vars.u2(X, A1, U1, A2, U2);
    std::cout << "u2=" << _u12 << std::endl;

    vec _u123 = vars.u3(X, A1, U1, A2, U2, A3, U3);
    std::cout << "u3=" << _u123 << std::endl;

    const int rep= 100000;
    boost::timer time;
    for(int i=0; i< rep; i++) {
	vec _u1234(vars.u4(X, A1, U1, A2, U2, A3, U3, A4, U4));
	//       std::cout << "u4=" << _u1234 << std::endl;
    }
    std::cout << "Compute time = " << 1000000.*time.elapsed() / rep << "us" << std::endl;
}
