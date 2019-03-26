// Filename: matrix_free_1.cpp (part of MTL4)

#include <iostream>
#include <cassert>
#include <boost/numeric/mtl/mtl.hpp>

struct poisson2D_dirichlet
{
    poisson2D_dirichlet(int m, int n) : m(m), n(n) {}

    template <typename Vector>
    Vector operator*(const Vector& v) const
    {
	mtl::vampir_trace<9901> tracer;
	assert(int(size(v)) == m * n);
	Vector w(m * n);
	
	for (int i= 0; i < m; i++)
	    for (int j= 0; j < n; j++) {
		int k= i * n + j; // offset
		w[k]= 4 * v[k];
		if (i > 0) w[k]-= v[k-n];   // upper neighbor
		if (i < m-1) w[k]-= v[k+n]; // lower neighbor
		if (j > 0) w[k]-= v[k-1];   // left neighbor
		if (j < n-1) w[k]-= v[k+1]; // right neighbor
	    }
	return w;
    }
    int m, n;
};

namespace mtl { namespace ashape {
    template <> struct ashape_aux<poisson2D_dirichlet> 
    {	typedef nonscal type;    };
}}

int main(int, char**)
{
    using namespace std;

    mtl::compressed2D<double> B;
    laplacian_setup(B, 1000, 1000);

    mtl::dense_vector<double> v(1000000), w(1000000);
    iota(v);

    poisson2D_dirichlet A(1000, 1000);
    // cout << "A * v is " << A * v << endl;
    w= A * v;
    
    w= B * v;


    return 0;
}
