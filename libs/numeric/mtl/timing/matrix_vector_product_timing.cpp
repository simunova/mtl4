#include <iostream>
#include <stdio.h>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/timer.hpp>

/*
  First run: 420ns dynamic types 
             369ns static types (r6809) 
             302ns static types with unrolled by hand (r6813)

	     1.1ns for v= A * u in r8527 -> switch to v+= 
	3x3: 8.9ns in r8527 with gen_mult v+= A * u
	     29ns dynamic
	     8.9ns mat_cvec
	     23.9ns mat_cvec dynamic
        5x5: 16ns 
	     73.3ns dynamic
	     16ns mat_cvec
	     71.2ns mat_cvec dynamic
*/

#define STATIC_TYPES

using namespace mtl;

#ifdef STATIC_TYPES
   typedef dense_vector<double, parameters<tag::col_major, fixed::dimension<3>, true> > vec;
   typedef mat::parameters<tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<3, 3>, true> mat_para;
   typedef dense_vector<double, parameters<tag::col_major, fixed::dimension<5>, true> > vec5;
   typedef mat::parameters<tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<5, 5>, true> mat_para5;
   #define VEC_ARG
#else
   typedef dense_vector<double> vec;
   typedef dense_vector<double> vec5;
   typedef mat::parameters<> mat_para;
   typedef mat::parameters<> mat_para5;
   #define VEC_ARG (3)
#endif
   typedef dense2D<double, mat_para> mat;
   typedef dense2D<double, mat_para5> mat5;


int main()
{
  {
      vec u(3), v(3, 0.0);
    u= 3., 4, 6; 

    mat A(3, 3);
    A=  2, 4, 8,
	8, 9, 1,
        4, 2, 1;

    const int rep= 100000000;
    boost::timer time;
    asm("#before loop");
    // for(int i= 0; i < rep; i++) 
	v+= A * u;
    
    asm("#after loop");
    std::cout << "Compute time = " << 1000000000.*time.elapsed() / rep << "ns" << std::endl;
    std::cout << "v is " << v << std::endl;
  }
  

  if (1) {
      vec5 u(5), v(5, 0.0);
    u= 3., 4, 6, 7, 9; 

    mat5 A(5, 5);
    A=  2, 4, 8, 0, 8,
	8, 9, 1, 2, 1,
	8, 9, 1, 2, 1,
	8, 9, 1, 2, 1,
        4, 2, 1, 4, 9;
    std::cout << "Nach Initialisierung\n";

    const int rep= 100000000;
    boost::timer time;
    for(int i= 0; i < rep; i++) {
	v+= A * u;
    }
    std::cout << "Compute time = " << 1000000000.*time.elapsed() / rep << "ns" << std::endl;
    std::cout << "v is " << v << std::endl;
  }
  return 0;
}
