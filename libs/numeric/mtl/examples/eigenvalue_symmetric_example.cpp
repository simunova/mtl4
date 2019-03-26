// Filename: eigenvalue_example.cpp (part of MTL4)

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    dense_vector<double>                    eig;

    double array[][4]= {{1,  1,   1,  0},
                        {1, -1,  -2,  0},
                        {1, -2,   1,  0},
                        {0,  0,   0, 10}};
    dense2D<double> A(array);
    std::cout << "A=\n" << A << "\n";

    eig= eigenvalue_symmetric(A,22);
    std::cout<<"eigenvalues  ="<< eig <<"\n";
    
    eig= 0;
    eig= qr_sym_imp(A);
    std::cout<<"eigenvalues  ="<< eig <<"\n";
    
    eig= 0;
    eig= qr_algo(A, 5);  // only 5 qr iterations (Q-R-changes)
    std::cout<<"eigenvalues  ="<< eig <<"\n";
 
    return 0;
}
