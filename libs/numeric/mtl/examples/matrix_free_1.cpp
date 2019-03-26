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

    mtl::dense_vector<double> v(20);
    iota(v);

    poisson2D_dirichlet A(4, 5);
    cout << "A * v is " << A * v << endl;

    return 0;
}
