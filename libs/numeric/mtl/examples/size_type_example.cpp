// Filename: size_type_example.cpp (part of MTL4)

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main()
{
    using namespace mtl;

    typedef mat::parameters<row_major, mtl::index::c_index, non_fixed::dimensions, false, unsigned> para;
    compressed2D<double, para>   A;
    laplacian_setup(A, 2, 3);

    std::cout << "A is\n" << A << "\nsize of index is " << sizeof(A.ref_minor()[0]) << '\n';
    return 0;
}
