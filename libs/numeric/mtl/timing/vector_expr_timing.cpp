#include <iostream>
#include <unistd>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/timer.hpp>

/*
  First run: 87ns + 74ns dynamic types 
             67ns + 42ns static types (r6809) (CET needs 52ns if run alone???)
             71ns + 46ns static types with unrolled by hand (r6810)
*/

#define STATIC_TYPES

#ifdef STATIC_TYPES
   typedef mtl::dense_vector<double, mtl::parameters<mtl::tag::col_major, mtl::fixed::dimension<3>, true> > vec;
#else
   typedef mtl::dense_vector<double> vec;
#endif



using namespace std;

int main()
{
    
    vec u(3), v(3), w(3), x(3);

    u= 3., 4, 6; v= 7, 9, 3; w= 7, 2, 4;

    const int rep= 10000000;
    boost::timer time;
    for(int i= 0; i < rep; i++) {
	x= dot(v, u) * w + 4.0 * v + 2 * w;
    }
    std::cout << "Compute time (CET) = " << 1000000000.*time.elapsed() / rep << "ns" << std::endl;

    time.restart();
    for(int i= 0; i < rep; i++) {
	x= dot(v, u) * (w+= 4.0 * v + 2 * w);
    }
    std::cout << "Compute time (RET) = " << 1000000000.*time.elapsed() / rep << "ns" << std::endl;

    return 0;
}
