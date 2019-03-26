#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    double           array[][3]= {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};
    dense2D<double>  A(array), A2, A3;

    // Creating a permutation matrix from a vector (or an array respectively)
    int indices[]= {1, 2, 0};
    mat::traits::permutation<>::type P= mat::permutation(indices);
    std::cout << "\nP =\n" << P;    

    // Permutating rows
    A2= P * A;
    std::cout << "\nP * A =\n" << A2;
    
    // Permutating columns
    A3= A2 * trans(P);
    std::cout << "\nA2 * trans(P) =\n" << A3;

    dense_vector<double> v(array[2]), w(P * v);
    std::cout << "\nP * v =\n" << w << "\n";
    
    return 0;
}
