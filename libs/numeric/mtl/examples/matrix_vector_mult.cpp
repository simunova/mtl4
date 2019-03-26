#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl; using namespace mtl::mat;
    
    const unsigned                xd= 2, yd= 5, n= xd * yd;
    dense2D<double>               A(n, n);
    compressed2D<double>          B(n, n);
    hessian_setup(A, 3.0); laplacian_setup(B, xd, yd); 

    typedef std::complex<double>  cdouble;
    dense_vector<cdouble>         v(n), w(n);
    for (unsigned i= 0; i < size(v); i++)
	v[i]= cdouble(i+1, n-i), w[i]= cdouble(i+n);

    v+= A * w;
    w= B * v;

    std::cout << "v is " << v << "\n";
    std::cout << "w is " << w << "\n";

    return 0;
}
