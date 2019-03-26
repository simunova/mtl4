// Filename: move_example.cpp (part of MTL4)

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

inline mtl::dense2D<double> make_identity(int m, int n)
{
    mtl::dense2D<double>  A(m, n);
    A= 1.0;
    return A;
}

int main()
{
    mtl::dense2D<double>   A;
    A= make_identity(3, 4);

    std::cout << "A is\n" << A;
    return 0;
}
