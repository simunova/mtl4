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

#include <string>
#include <iostream>

#include <boost/numeric/mtl/mtl.hpp>


using namespace std;


template <typename Matrix, typename Tag>
void two_d_iteration(char const* outer, const Matrix& matrix, Tag)
{
    namespace traits = mtl::traits;

    typename traits::row<Matrix>::type                                 row(matrix); 
    typename traits::col<Matrix>::type                                 col(matrix); 
    typename traits::const_value<Matrix>::type                         value(matrix); 
    typedef typename traits::range_generator<Tag, Matrix>::type        cursor_type;

    cout << outer << '\n';
    for (cursor_type cursor = mtl::begin<Tag>(matrix), cend = mtl::end<Tag>(matrix); cursor != cend; ++cursor) {
	typedef mtl::tag::nz     inner_tag;
	cout << "---\n";
	typedef typename traits::range_generator<inner_tag, cursor_type>::type icursor_type;
	for (icursor_type icursor = mtl::begin<inner_tag>(cursor), icend = mtl::end<inner_tag>(cursor); icursor != icend; ++icursor)
	    cout << "matrix[" << row(*icursor) << "][" << col(*icursor) << "] = " << value(*icursor) << '\n';
    }
} 


template <typename Matrix, typename Tag>
void two_d_iteration(char const* name, const Matrix&, Tag, mtl::complexity_classes::infinite)
{
    cout << name << ": Tag has no implementation\n";
}

 
template <typename Matrix, typename Value>
void test(string name, const Matrix& A, Value check)
{
    cout << name << '\n';
    cout << "num_rows(A) is " << num_rows(A) << '\n';
    cout << "num_cols(A) is " << num_cols(A) << '\n';
    cout << "A[0][0] is " << A[0][0] << '\n';
    cout << name << ", A is\n" << A << '\n';
    
    MTL_THROW_IF(A[0][1] != check, mtl::runtime_error("Wrong value in A[0][1]!"));

    two_d_iteration("Row-wise", A, glas::tag::row());
    two_d_iteration("Column-wise", A, mtl::tag::col());
    two_d_iteration("On Major", A, mtl::tag::major());


    mtl::transposed_view<const Matrix> At(A);
    cout << "\n===\nA^T is\n" << At << '\n';
    MTL_THROW_IF(At[1][0] != check, mtl::runtime_error("Wrong value in At[1][0]!"));

    two_d_iteration("Transposed row-wise", At, mtl::tag::row());
    two_d_iteration("Transposed Column-wise", At, mtl::tag::col());
    two_d_iteration("Transposed On Major", At, mtl::tag::major());

    mtl::dense2D<double> B(3.0 * A);
    cout << "3*A is\n" << B << "\n\n";
}

int main(int, char**)
{
    typedef mtl::dense_vector<double> vt;
    typedef mtl::dense_vector<double, mtl::vec::parameters<mtl::row_major> >  vrt;
    vt  v(2), w(3);
    vrt z(3);
    v= 1, 2; w= 2, 3, 4; z= trans(w);

#if 0
    using namespace mtl; 
    cout << "ashape<v> is " << typeid(ashape::ashape<vt>::type).name() << "\n";
    cout << "ashape<z> is " << typeid(ashape::ashape<vrt>::type).name() << "\n";
    cout << "ashape<v*z> is " << typeid(ashape::mult_op<ashape::ashape<vt>::type, ashape::ashape<vrt>::type>::type).name() << "\n";
    cout << "mult_result<v*z> is " << typeid(mtl::traits::vec_mult_result<vt, vrt>::type).name() << "\n";
    cout << "v * z is " << typeid(v * z).name() << "\n";
    return 0;
#endif

    test("ones(2, 3)", mtl::ones(2, 3), 1);
    test("ones<float>(2, 3)", mtl::ones<float>(2, 3), 1.f);
    test("v * trans(w)", mtl::mat::outer_product_matrix<vt, vt>(v, w), 3.0);
    test("v * z", v * z, 3.0);
    test("v * trans(w)", v * trans(w), 3.0);
    test("hilbert_matrix(2, 3)", mtl::mat::hilbert_matrix<>(2, 3), 0.5);
    
    return 0;
}
 
