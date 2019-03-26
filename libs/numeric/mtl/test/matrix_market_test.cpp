// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#include <iostream>
#include <cmath>
#include <string>
#include <boost/type_traits/is_complex.hpp>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/recursion/matrix_recursator.hpp>
 
using namespace std; 
using namespace mtl;

std::string program_dir; // Ugly global variable !!!

template <typename Matrix>
void inline test_file(Matrix& A, const char* file_name, const char* comment)
{
    cout << "Filename is " << mtl::io::join(program_dir, file_name) << "\n";
    mtl::io::matrix_market_istream ms(mtl::io::join(program_dir, file_name));
    ms >> A;
    std::cout << "Read from " << file_name << " (" << comment << ") is " 
	      << num_rows(A) << "x" << num_cols(A) << "\n";

    if (num_rows(A) > 9 && num_cols(A) > 9) {
	int reordering[]= {0, 1, 2, 3, 4, 5, 6, 7, 8};
	mtl::mat::traits::reorder<>::type  R= mtl::mat::reorder(reordering, num_cols(A)),
	                                      R2= mtl::mat::reorder(reordering, num_rows(A));
	Matrix B0(R * A), B(B0 * trans(R2));
	std::cout << "A[0:9][0:9] is:\n" << B;
    } else
	std::cout << "A is:\n" << A;
}



template <typename Matrix>
void inline test(Matrix& A, const char* name)
{
    std::cout << "\n" << name << "\n";

    typedef typename mtl::Collection<Matrix>::value_type vt;

    test_file(A, "matrix_market/jgl009.mtx", "general pattern");
    if (boost::is_complex<vt>::value)
	test_file(A, "matrix_market/mhd1280b.mtx", "Hermitian"); 
    // test_file(A, "matrix_market/plskz362.mtx", "Skew-symmetric"); // has only 0s in A[:9][:9]
    test_file(A, "matrix_market/bcsstk01.mtx", "Real symmetric");

    Matrix B(mtl::io::matrix_market(mtl::io::join(program_dir, "matrix_market/jgl009.mtx"))), C;
    std::cout << "Matrix market file read in constructor:\n" << B;

    C= mtl::io::matrix_market(mtl::io::join(program_dir, "matrix_market/jgl009.mtx"));
    std::cout << "Matrix market file assigned:\n" << B;
 }

template <typename Matrix>
void inline failure_test(Matrix& A)
{
    try {
	A= mtl::io::matrix_market("File_not_exist_test.mtx");
    } catch (const mtl::file_not_found& e) {
	std::cerr << "Successfully caught exception for inexistant file. Error message is:\n" << e.what();
	return;
    }
    throw "No exception thrown for inexistant file.";
}

template< typename Matrix >
void read_test(Matrix& A, const char* name )
{
	typedef typename Collection<Matrix>::size_type my_size;
	dense2D< double > data(2,4);
	unsigned v(1);
	for(unsigned r(0); r < 2; ++r)
		for(unsigned c(0); c < 4; ++c,++v)
			data(r,c)=v;
	mtl::io::matrix_market_istream in("matrix_market/dense_read.mtx");
	in >> A;
	dense2D< double > diff(A-data);
	pair< my_size, my_size > pos(max_abs_pos(diff));
	if(diff(pos.first,pos.second) > 1e-14) {
		std::stringstream ss;
		ss << "could not read the dense matrix with type " << name << std::endl;
		ss << "wanted:\t" << data << "\n";
		ss << "got:\t" << A << "\n";
		throw std::runtime_error(ss.str());
	}
}

int main(int, char* argv[])
{
    using namespace mtl;

    compressed2D<double>                             cdc;
    compressed2D<std::complex<double> >              ccc;
    dense2D<double>                                  dc;
    dense2D<double, mat::parameters<col_major> >  dcc;
    dense2D<float>                                   fc;
    morton_dense<double,  morton_mask>               mdc;
    morton_dense<double, doppled_32_col_mask>        mcc;

    program_dir= mtl::io::directory_name(argv[0]);

    test(cdc, "compressed2D");
    test(ccc, "compressed2D complex");
    test(dc, "dense2D");
    test(dcc, "dense2D col-major");
    test(mdc, "pure Morton");
    test(mcc, "Hybrid col-major");

    failure_test(cdc);
 
    read_test(dc, "dense2D");
    //read_test(ccc, "compressed2D complex");
    read_test(dcc, "dense2D col-major");
    read_test(mdc, "pure Morton");
    read_test(mcc, "Hybrid col-major");

    return 0;
}
