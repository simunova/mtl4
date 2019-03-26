// Filename: size_type_example2.cpp (part of MTL4)

#include <iostream>
#include <boost/cstdint.hpp>
#include <boost/numeric/mtl/mtl.hpp>

int main()
{
    using namespace mtl;

    typedef mat::parameters<row_major, mtl::index::c_index, non_fixed::dimensions, false, boost::uint_least32_t> para;
    compressed2D<double, para>   A;
    laplacian_setup(A, 2, 3);
    
    std::cout << "A is\n" << A << "\nsize of index is " << sizeof(A.ref_minor()[0]) << '\n';
    return 0;
}
