#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    double           array[][3]= {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};
    dense2D<double>  A(array), B2, B3;

    // Creating a reordering matrix from a vector (or an array respectively)
    int indices[]= {2, 1, 1, 2};
    mat::traits::reorder<>::type R= mat::reorder(indices);
    std::cout << "\nR =\n" << R;    

    // Reorder rows
    B2= R * A;
    std::cout << "\nR * A =\n" << B2;
    
    // Reorder columns
    B3= B2 * trans(R);
    std::cout << "\nB2 * trans(R) =\n" << B3;
    
    dense_vector<double> v(array[2]), w(R * v);
    std::cout << "\nR * v =\n" << w << "\n";
    
    return 0;
}
