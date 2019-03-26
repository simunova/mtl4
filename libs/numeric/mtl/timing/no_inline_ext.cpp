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
#include <vector>
#include <boost/timer.hpp>


#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>
#include <boost/numeric/mtl/operation/matrix_mult.hpp>
#include <boost/numeric/mtl/matrix/hessian_setup.hpp>
#include <boost/numeric/mtl/operation/assign_mode.hpp>


using namespace mtl;
using namespace mtl::recursion; 
using namespace std;  


typedef gen_tiling_44_dense_mat_mat_mult_t<assign::plus_sum>  tiling_44_base_mult_t;

#if 0
    // Bitmasks: 
    const unsigned long 
	doppled_64_row_mask= generate_mask<true, 6, row_major, 0>::value,
	doppled_64_col_mask= generate_mask<true, 6, col_major, 0>::value;

void hybrid_ext_mult_44(const morton_dense<double,  doppled_64_row_mask>& a, 
			const morton_dense<double,  doppled_64_col_mask>& b,
			morton_dense<double,  doppled_64_row_mask>& c)
{
  tiling_44_base_mult_t()(a, b, c);
}
#endif

void dense_ext_mult_44(const dense2D<double>& a,
		       const dense2D<double, mat::parameters<col_major> >& b,
		       dense2D<double>& c)
{
  tiling_44_base_mult_t()(a, b, c);
}
 
