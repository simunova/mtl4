#include <iostream>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/timer.hpp>

#ifdef __AVX__ 
#include <immintrin.h>
#include <cstdlib>
#endif

const unsigned rep= 100000000; 

template <unsigned BSize, typename Vector>
void run(Vector& u, const Vector& v, Vector& w, Vector&)
{
    using mtl::unroll;

    boost::timer t1; 
    for (unsigned j= 0; j < rep; j++)
	unroll<BSize>(u)= v + w;
    std::cout << "Compute time unroll<" << BSize << ">(u)= v + v is " 
	      << 1000000.0 * t1.elapsed() / double(rep) << " us = "
	      << double(rep) * size(u) / t1.elapsed() / 1e9  << " GFlops.\n";
    // std::cout << "u is " << u << '\n';

}

struct add
{
    template <typename T>
    T operator()(const T& x, const T& y)
    { return x + y; }
};

struct as
{
    template <typename T>
    T& operator()(T& x, const T& y)
    { return x= y; }
};

template <typename T>
void f(T& a)
{
    T& x= a;
    x= "fasdf";
}

int main(int argc, char** argv)
{
    using namespace mtl;
    const unsigned cs= 50;
    unsigned s= cs;
    if (argc > 1) s= atoi(argv[1]);
    
#if 0
    double a[cs] __attribute__ ((aligned (__BIGGEST_ALIGNMENT__))), 
           b[cs] __attribute__ ((aligned (__BIGGEST_ALIGNMENT__))), 
	   c[cs] __attribute__ ((aligned (__BIGGEST_ALIGNMENT__)));
#else
    double *a= new double[s], 
           *b= new double[s], 
	   *c= new double[s];
#endif

#if 0
    typedef dense_vector<double, parameters<tag::col_major, fixed::dimension<cs>, true > > vec;
    vec u, v, w, x;
#else
    typedef mtl::dense_vector<double> vec;
    double *ua, *va, *wa, *xa;
    posix_memalign( (void**)&ua, 32, 8 * s);
    posix_memalign( (void**)&va, 32, 8 * s);
    posix_memalign( (void**)&wa, 32, 8 * s);
    posix_memalign( (void**)&xa, 32, 8 * s);
    vec u(s, ua), v(s, va), w(s, wa), x(s, xa);
#endif   

    for (unsigned i= 0; i < s; i++) { 
	a[i]= v[i]= double(i);
	b[i]= w[i]= double(2*i + 15);
    }
    boost::timer t; 
    asm("#before loop");
#ifdef __AVX__ 
    __m256d *up= reinterpret_cast<__m256d*>(&u[0]);
    __m256d *vp= reinterpret_cast<__m256d*>(&v[0]);
    __m256d *wp= reinterpret_cast<__m256d*>(&w[0]);
    const unsigned cs4= s/4;
    for (unsigned j= 0; j < rep; j++) {
	for (unsigned i= 0; i < cs4; i++) 
	    up[i]= vp[i] + wp[i];
	for (unsigned i= 4 * cs4; i < s; i++) 
	    u[i]= v[i] + w[i];
    }
#else
    // #warning "No AVX!"
    for (unsigned j= 0; j < rep; j++) 
	 u= v + w;
#endif
    asm("#after loop");
    std::cout << "Compute time u= v + v is " << 1000000.0 * t.elapsed() / double(rep) << " us = "
	      << double(rep) * size(u) / t.elapsed() / 1e9  << " GFlops.\n";


#if 0
    //double *&ar= a;
    double (&ar)[cs] __attribute__ ((aligned (__BIGGEST_ALIGNMENT__)))= reinterpret_cast<double (&)[cs]>(v[0]); // perverser cast auf array ref
    double (&br)[cs] __attribute__ ((aligned (__BIGGEST_ALIGNMENT__)))= reinterpret_cast<double (&)[cs]>(w[0]);
    double (&cr)[cs] __attribute__ ((aligned (__BIGGEST_ALIGNMENT__)))= reinterpret_cast<double (&)[cs]>(u[0]);
#endif
    // f(a);

    add adder;
    as  assigner;
    t.restart();
    asm("#before c-loop");
    for (unsigned j= 0; j < rep; j++)
	for (unsigned i= 0; i < cs; i++) {
	    // assigner(cr[i], adder(ar[i], br[i]));
	    assigner(c[i], adder(a[i], b[i]));
	}
    asm("#after c-loop");
    std::cout << "Compute time u= v + v is " << 1000000.0 * t.elapsed() / double(rep) << " us = "
	      << double(rep) * s / t.elapsed() / 1e9  << " GFlops.\n";

#if 0

    run<1>(u, v, w, x);
    run<2>(u, v, w, x);
    run<4>(u, v, w, x);
    run<6>(u, v, w, x);
    run<8>(u, v, w, x);
#endif

    std::cout << "c[0] = " << c[0] << ", u[0] = " << u[0] << '\n';
    c[0]= a[0] + b[0];
    return 0 ;
}

