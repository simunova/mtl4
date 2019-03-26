#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;
    
    typedef std::complex<double>      cdouble;
    const unsigned                    xd= 2, yd= 5, n= xd * yd;
    compressed2D<cdouble>             A(n, n);
    mat::laplacian_setup(A, xd, yd); 

    // Fill imaginary part of the matrix
    A*= cdouble(1, -1);
    std::cout << "A is\n" << with_format(A, 7, 1) << "\n";

    std::cout << "trace(A) is " << trace(A) << "\n\n";
    std::cout << "conj(A) is\n" << with_format(mtl::mat::conj(A), 7, 1) << "\n"; // ADL issue on g++ 4.4
    std::cout << "trans(A) is\n" << with_format(trans(A), 7, 1) << "\n";
    std::cout << "hermitian(A) is\n" << with_format(hermitian(A), 7, 1) << "\n";

    return 0;
}
