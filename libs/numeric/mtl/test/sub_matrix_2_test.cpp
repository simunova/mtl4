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
#include <boost/numeric/mtl/recursion/matrix_recursator.hpp>


using namespace std;  


template <typename Matrix>
void test(Matrix& A, const char* name)
{
    using mtl::irange;
    A= 0.0;
    A[0][0]= 1.0; 
    hessian_setup(A, 1.0);

    A[irange(0, 8)][irange(0, 8)];
    mtl::mat::recursator<Matrix> rec(A);

    std::cout << "\n" << name << "\n";
    std::cout << "A:\n" << A << '\n';    
    std::cout << "A[irange(0, 8)][irange(0, 8)]:\n" << A[irange(0, 8)][irange(0, 8)] << '\n';    
    std::cout << "*rec:\n" << *rec << '\n';    

    mtl::mat::recursator<Matrix> nw= north_west(rec);
    std::cout << "north_west:\n" << *nw << '\n';    

    std::cout << "north_west of north_west:\n" << *north_west(nw) << '\n';
    MTL_THROW_IF((*north_west(nw))[0][0] != 0.0, mtl::runtime_error("(*north_west(nw))[0][0] != 0.0"));
    (*north_west(nw))[0][0]= 2.0;

    std::cout << "south_east of north_west:\n" << *south_east(nw) << '\n';
    MTL_THROW_IF((*south_east(nw))[0][0] != 4.0, mtl::runtime_error("(*south_east(nw))[0][0] != 4.0"));

    std::cout << "north_west of north_west:\n" << *north_west(nw) << '\n';
    MTL_THROW_IF((*north_west(nw))[0][0] != 2.0, mtl::runtime_error("(*north_west(nw))[0][0] != 2.0"));

    std::cout << "south_east of north_west:\n" << *south_east(nw) << '\n';
    MTL_THROW_IF((*south_east(nw))[0][0] != 4.0, mtl::runtime_error("(*south_east(nw))[0][0] != 4.0"));

    std::cout << "nw.first_address() == " << nw.first_address() 
	      << ", &(*nw)[0][0] == " << &(*nw)[0][0] << '\n';
    MTL_THROW_IF(nw.first_address() != &(*nw)[0][0], mtl::runtime_error("Inconsistency in address calculation"));
}


int main(int, char**)
{
    using namespace mtl;
    const unsigned size= 5; 

    dense2D<double> dc(size, size-2);
    dense2D<double, mat::parameters<col_major> >  dcc(size, size-2);
    dense2D<float>                                   fc(size, size-2);
    morton_dense<double,  morton_mask>               mdc(size, size-2);
    morton_dense<double, doppled_32_col_mask>        mcc(size, size-2);

    test(dc, "dense2D");
    test(dcc, "dense2D col-major");
    test(mdc, "pure Morton");
    test(mcc, "Hybrid col-major");

    return 0;
}
