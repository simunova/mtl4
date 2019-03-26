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
#include <boost/numeric/mtl/mtl.hpp>

template <typename ArgMatrix>
void inline f(const ArgMatrix&) {}


template <typename ResMatrix, typename ArgMatrix>
void test(const ResMatrix& , const ArgMatrix& B)
{
    ResMatrix C(trans(B));

    std::cout << "trans(B) is \n" << C;
    std::cout << "trans(B) is \n" << trans(B);
    std::cout << "trace(B) is " << trace(B) << "\n";

    // ResMatrix D(A * trans(B));

#if 0 // Can't convert so far
    ResMatrix D= B;
    std::cout << "D(B) is \n" << D;
    
    f<ResMatrix>(B);
#endif
}


int main(int , char**)
{
    using namespace mtl;
    typedef mat::parameters<tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<2, 2>, true> fmat_para;

    float ma[2][2]= {{2., 3.}, {4., 5.}};
    
    dense2D<float>                   A_dyn(ma);
    dense2D<float, fmat_para>        A_stat(ma);

    test(A_dyn, A_dyn);
    test(A_dyn, A_stat);
    test(A_stat, A_stat);
    test(A_stat, A_dyn);

    typedef mtl::vec::fixed::dimension< 2 > fsize;
    mtl::dense_vector<float, mtl::vec::parameters<mtl::col_major, fsize, true> >     rf;
    rf = 1.0;
    mtl::dense_vector<float, mtl::vec::parameters<mtl::col_major, fsize, true> >     rf2;
    const static double swap = 1.0;
    swap*(rf + A_stat*rf);
    return 0;
}

