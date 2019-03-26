#include <gmpxx.h>
#include <boost/test/minimal.hpp>
#include <boost/numeric/mtl/mtl.hpp>

// Check only if malloc error happens as reported by Hui Li

int test_main(int argc, char** argv) 
{
    mtl::dense2D<mpq_class> A(2, 2);
    A(0, 0) = mpq_class(0);
    A(0, 1) = 0; // problem
    mpq_class x;
    x = 0;
    return 0;
}



// g++ -g -DMTL_ASSERT_FOR_THROW -I$MTL_BOOST_ROOT -I/home/pgottsch/Download/gmp-4.2.2 -I/home/pgottsch/projects/boost/boost_1_33_1 -o gmp_test gmp_test.cpp -L/home/pgottsch/Download/gmp-4.2.2 -lgmp

