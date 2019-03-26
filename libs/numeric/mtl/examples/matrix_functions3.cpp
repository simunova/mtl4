#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl; using mtl::iall;

    typedef std::complex<double>      cdouble;
    const unsigned                    xd= 2, yd= 5, n= xd * yd;
    dense2D<cdouble>                  A(n, n);
    mat::laplacian_setup(A, xd, yd); 

    // Fill imaginary part of the matrix
    A*= cdouble(1, -1);
    std::cout << "A is\n" << with_format(A, 7, 1) << "\n";

    std::cout << "sub_matrix(A, 2, 4, 1, 7) is\n" 
	      << with_format(sub_matrix(A, 2, 4, 1, 7), 7, 1) << "\n";

    //col-vector from matrix
    dense_vector<cdouble>   v_c(A[iall][0]);

    std::cout << "col-vector v_c is\n" << v_c << "\n";

    //row-vector from matrix
    dense_vector<cdouble, mtl::vec::parameters<tag::row_major> > v_r(A[0][iall]);

    std::cout << "row-vector v_r is\n" << v_r << "\n";

    //row-vector in matrix
    RowInMatrix<dense2D<cdouble> >::type v_r2(A[0][iall]);

    std::cout << "row-vector v_r2 is\n" << v_r2 << "\n";

    //submatrix from matrix per begin and end of row and column
    dense2D<cdouble> B= sub_matrix(A, 2, 4, 1, 7);
    B[1][2]= 88;

    std::cout << "B is\n" << B << "\n";

    //submatrix from matrix per irange
    using mtl::irange;
    irange row(2, 4), col(1, 7);
    dense2D<cdouble> B1= A[row][col];

    std::cout << "B1 is\n" << B1 << "\n";

    //scalar from matrix
    cdouble C= A[1][1];

    std::cout << "C is\n" << C << "\n";

    return 0;
}
