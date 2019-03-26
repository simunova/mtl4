#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl; using namespace mtl::mat;
    
    const unsigned                xd= 2, yd= 5, n= xd * yd;
    dense2D<double>               A(n, n);
    laplacian_setup(A, xd, yd); 
    dense_vector<double>          v(n), w(n, 7.0);

    // Scale A with 4 and multiply the scaled view with w
    v= 4 * A * w;
    std::cout << "v is " << v << "\n";

    // Scale w with 4 and multiply the scaled view with A
    v= A * (4 * w);
    std::cout << "v is " << v << "\n";

    // Scale both with 2 before multiplying
    v= 2 * A * (2 * w);
    std::cout << "v is " << v << "\n";

    // Scale v after the MVP
    v= A * w;
    v*= 4;
    std::cout << "v is " << v << "\n";

    return 0;
}
