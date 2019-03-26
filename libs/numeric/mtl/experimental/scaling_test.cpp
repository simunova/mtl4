#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;

template <int Scale> struct scaling 
{    template <typename T> T operator*(T x) const { return Scale * x; }  };

template <> struct scaling<-1> 
{    template <typename T> T operator*(T x) const { return -x; }   };

template <> struct scaling<0>
{    template <typename T> T operator*(T x) const { return T(0); } };

template <> struct scaling<1>
{    template <typename T> T operator*(T x) const { return x; }    };


typedef mtl::dense_vector<int,   mtl::parameters<mtl::col_major, mtl::fixed::dimension<2>, true> > fveci2_type;
typedef mtl::dense_vector<float, mtl::parameters<mtl::col_major, mtl::fixed::dimension<2>, true> > fvecf2_type;

template <typename S, typename V> struct scaled_product2;

template <int S0, int S1>
struct scaling_vector2 : public fveci2_type
{
    typedef scaling<S0>             s0;
    typedef scaling<S1>             s1;
    typedef scaling_vector2<S0, S1> self;

    template <typename Vector>
    scaled_product2<self, Vector> operator*(const Vector& v)
    {
	return scaled_product2<self, Vector>(v);
    }
};

template <typename S, typename V>
struct scaled_product2 : public fvecf2_type
{
    explicit scaled_product2(const V& v) : v(v) {}

    typedef typename mtl::Collection<V>::value_type value_type;
    typedef typename mtl::Collection<V>::size_type  size_type;

    value_type operator[](size_type i) 
    {
	return i == 0 ? typename S::s0() * v[0] : typename S::s1() * v[1];
    }
    
    // template <int I> value_type at() { return S::at<I>::type() * v[I]; }

    const V& v;
};

struct integrator
{
    /* const */ static fvecf2_type s0;
};

float s0a[2]= {1.0f, -1.0f};
fvecf2_type integrator::s0(s0a); // (*this)[0]= };

int main(int argc, char** argv) 
{
    scaling<-1> c;
    scaling<1>  d;

    std::cout << "c * 3.0 + d * 4.5f = " << c * 3.0 + d * 4.5f << '\n';

    fvecf2_type v, w, z;
    v= 3.0, 4.5;
    w= 2.0, 7.5;
    std::cout << "v = " << v << ", w = " << w << '\n';
#if 0
    scaling_vector2<-1, 1> s0;
    scaling_vector2< 1, 0> s1;

    z= s0 * v; // + s1 * w;
    std::cout << "z = " << z << '\n';
#endif
    return 0;
}
