// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschrÃ¤nkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>


using namespace std;

int main(int , char** ) 
{
  
    typedef mtl::mat::parameters<mtl::row_major,
	mtl::index::c_index,
	mtl::non_fixed::dimensions,
	false, int> para;
    typedef mtl::compressed2D<double,para> matrix_type;
    matrix_type Ks;
    //
    mtl::io::matrix_market_istream ik("../../../../../branches/data/matrix_market/outK.mm");
    //
    ik >> Ks;
    //
    ik.close();
    std::cout << num_rows(Ks) << std::endl;
    std::cout << num_cols(Ks) << std::endl;
    std::cout << Ks.nnz() << std::endl;
    for(int k=0;k<10;k++)
        std::cout << Ks(k,k) << std::endl; 

    if (num_rows(Ks) < 20)
	std::cout << Ks;

    return 0;
}
