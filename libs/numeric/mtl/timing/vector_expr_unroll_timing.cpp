#include <iostream>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/timer.hpp>

/*
  First run: 87ns + 74ns dynamic types 
             67ns + 42ns static types (r6809) (CET needs 52ns if run alone???)
             71ns + 46ns static types with unrolled by hand (r6810)
*/

const unsigned rep= 1000000; 

template <unsigned BSize, typename Vector>
void run(Vector& u, const Vector& v, Vector& w, Vector& x)
{
    using mtl::unroll;

    boost::timer t1; 
    for (unsigned j= 0; j < rep; j++)
	unroll<BSize>(u)= v + v + w;
    std::cout << "Compute time unroll<" << BSize << ">(u)= v + v + w is " 
	      << 1000000.0 * t1.elapsed() / double(rep) << " 탎.\n";
    std::cout << "u is " << u << '\n';

    boost::timer t2; 
    for (unsigned j= 0; j < rep; j++)
	unroll<BSize>(x)= dot(v, u) * (w+= 4.0 * v + 2 * w);
    std::cout << "Compute time unroll<" << BSize << ">(x)= dot(v, u) * (w+= 4.0 * v + 2 * w) is " 
	      << 1000000.0 * t2.elapsed() / double(rep) << " 탎.\n";
    std::cout << "u is " << u << '\n';
}

int main(int argc, char** argv)  
{
    unsigned s= 1000;
    if (argc > 1) s= atoi(argv[1]);
    mtl::dense_vector<float> u(s), v(s), w(s), x(s);

    for (unsigned i= 0; i < s; i++) { 
	v[i]= float(i);
	w[i]= float(2*i + 15);
    }
    boost::timer t; 
    for (unsigned j= 0; j < rep; j++)
	u= v + v + w;
    std::cout << "Compute time u= v + v + w is " << 1000000.0 * t.elapsed() / double(rep) << " 탎.\n";
    std::cout << "u is " << u << '\n';

    boost::timer t2; 
    for (unsigned j= 0; j < rep; j++)
	x= dot(v, u) * (w+= 4.0 * v + 2 * w);
    std::cout << "Compute time x= dot(v, u) * (w+= 4.0 * v + 2 * w) is " << 1000000.0 * t2.elapsed() / double(rep) << " 탎.\n";
    std::cout << "u is " << u << '\n';

    run<1>(u, v, w, x);
    run<2>(u, v, w, x);
    run<4>(u, v, w, x);
    run<6>(u, v, w, x);
    run<8>(u, v, w, x);

    return 0 ;
}

