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
#include <boost/numeric/mtl/mtl.hpp>


template <typename ResMatrix, typename ArgMatrix>
void test(const ResMatrix&, const ArgMatrix& B)
{
    ResMatrix C(B * B);

    C+= trans(B) * B;
    C+= trans(B) * B * B;
#if 0
	std::cout << typeid(typename mtl::traits::category<mtl::mat::mat_mat_times_expr<ArgMatrix, ArgMatrix> >::type).name() << '\n';
	std::cout << typeid(typename mtl::traits::category<mtl::mat::rscaled_view<ArgMatrix, double> >::type).name() << '\n';
	char c; std::cin >> c; 
#endif
	C+= B * 3.5 * B * B;
    C+= trans(B) * 3.5 * B * B;

    C+= 3.5 * ArgMatrix(B * B);
    C= 3.5 * ArgMatrix(B * B);

    //C+= 3.5 * (B * B);
    //C= 3.5 * (B * B);
}


int main(int, char**)
{
    using namespace mtl;
    typedef mat::parameters<tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<2, 2>, true> fmat_para;

    float ma[2][2]= {{2., 3.}, {4., 5.}};
    
    dense2D<float>                   A_dyn(ma);
    dense2D<float, fmat_para>        A_stat(ma);

    test(A_dyn, A_dyn);
    //test(A_dyn, A_stat);

    return 0;
}

