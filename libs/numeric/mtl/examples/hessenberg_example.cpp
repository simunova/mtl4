#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;
    dense2D<double>       A(5, 5), H(5, 5);  A= 0.0;  
    
    A[0][1] = 3;  A[1][4] = 7; A[0][0] = 1; A[4][4] = 17;
    A[2][3] = -2; A[2][4] = 5; A[4][0] = 2; A[4][1] = 3;
    A[3][2] = 4;
    
    H= hessenberg(A);
    std::cout<< "Hessenberg=\n" << H << "\n";
    H= extract_householder_hessenberg(A);
    std::cout<< "extract_householder_hessenberg=\n" << H << "\n";
    H= extract_hessenberg(A);
    std::cout<< "extract_hessenberg=\n" << H << "\n";
    H= householder_hessenberg(A);
    std::cout<< "householder_hessenberg=\n" << H << "\n";
    H= hessenberg_factors(A);
    std::cout<< "hessenberg_factors=\n" << H << "\n";
    // H= hessenberg_q(A);
    // std::cout<< "hessenberg_q=\n" << H << "\n";

   return 0;
}
