#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    double           array[][3]= {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};
    dense2D<double>  A(array);

    std::cout << "\nupper(A) = \n" << upper(A);

    std::cout << "\nstrict_upper(A) = \n" << strict_upper(A);
        
    return 0;
}
