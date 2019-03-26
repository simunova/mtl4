#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;
int main(int, char**)
{
    mtl::dense2D<double>    A(3, 4), A_t(3, 4), S(3, 3), V(3, 4), D(4,4);
    A= 0;

    A[0][0]=1;  A[0][1]=1;  A[0][2]=1;  A[0][3]=4; 
    A[1][0]=1;  A[1][1]=2;  A[1][2]=2;  A[1][3]=3;
    A[2][0]=9;  A[2][1]=3;  A[2][2]=2;  A[2][3]=4;       
    std::cout<<"A=\n"<< A <<"\n";

    boost::tie(S, V, D)= svd(A, 1.e-10)= svd(A);  //second argument is optional (missmatch of upper R (A= Q*R))
    std::cout<<"Matrix  S=\n"<< S <<"\n";
    std::cout<<"Matrix  V=\n"<< V <<"\n";
    std::cout<<"Matrix  D=\n"<< D <<"\n";
    A_t= S*V*trans(D);
    std::cout<<"Matrix  A=S*V*D'=\n"<< A_t <<"\n";

    return 0;
}

