#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int , char**)
{
    using namespace mtl;
    typedef mtl::vec::parameters<tag::col_major, mtl::vec::fixed::dimension<2> > fvec_para;
    typedef mat::parameters<tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<2, 2> > fmat_para;

    dense2D<float, fmat_para>        A, B; // dimension not needed here
    dense_vector<float, fvec_para>   v, w; // here neither

    A= 2., 3.,
       4., 5.;
    v= 3., 4.;

    w= A * v; // Same syntax as dynamic size
    B= A * A; 

    std::cout << "A * v is " << w << "\n\n";
    std::cout << "A * A is\n" << B;

    return 0;
}
