// Software License for PMTL
// 
// Copyright (c) 2010 SimuNova UG, www.simunova.com.
// All rights reserved.
// Author: Peter Gottschling
// 
// This file is part of the Parallel Matrix Template Library
// 
// The details are regulated by the EULA at http://www.simunova.com/en/eula
//                             respectively http://www.simunova.com/de/agb.

#include <iostream>
#include <typeinfo>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/utility/flatcat.hpp>

using namespace std;  
using namespace mtl;

template <typename U>
inline void check(const U&) 
{
    io::tout << "correctly dispatched to " << typeid(U).name() << '\n';
}


template <typename T, typename U1, typename U2, typename U3>
void test(const T&, const U1&, const U2&, const U3&)
{
    io::tout << "Type " << typeid(T).name() << ":\n - trying to match <dense> ";
    check<U1>(traits::flatcat1<T, tag::dense>());
    io::tout << " - trying to match <multi_vector, dense> ";
    check<U2>(traits::flatcat2<T, tag::multi_vector, tag::dense>());
    io::tout << " - trying to match <multi_vector, dense, sparse> ";
    check<U3>(traits::flatcat3<T, tag::multi_vector, tag::dense, tag::sparse>());

    io::tout << "\n";
}



int main(int, char**) 
{
    typedef dense_vector<float>                        	     	      dvf;
    typedef dense_vector<std::complex<double> >        	     	      dvc;
    typedef dense_vector<float, mtl::vec::parameters<row_major> >  dvr;

    typedef dense2D<double>                                           dmr;
    typedef dense2D<double, mat::parameters<tag::col_major> >      dmc;
    typedef morton_dense<double, morton_mask>                         mmd;
    typedef compressed2D<double>                                      cmr;
    typedef compressed2D<double, mat::parameters<tag::col_major> > cmc;
    typedef multi_vector<dvf>                                         mvc;
    
    test(dvf(), tag::flat<tag::dense>()    , tag::flat<tag::dense>()	   , tag::flat<tag::dense>());
    test(dvc(), tag::flat<tag::dense>()	   , tag::flat<tag::dense>()	   , tag::flat<tag::dense>());
    test(dvr(), tag::flat<tag::dense>()	   , tag::flat<tag::dense>()	   , tag::flat<tag::dense>());

    test(dmr(), tag::flat<tag::dense>()	   , tag::flat<tag::dense>()	   , tag::flat<tag::dense>());
    test(dmc(), tag::flat<tag::dense>()	   , tag::flat<tag::dense>()	   , tag::flat<tag::dense>());
    test(mmd(), tag::flat<tag::dense>()	   , tag::flat<tag::dense>()	   , tag::flat<tag::dense>());
    test(cmr(), tag::universe()            , tag::universe()        	   , tag::flat<tag::sparse>());
    test(cmc(), tag::universe()            , tag::universe()        	   , tag::flat<tag::sparse>());
    test(mvc(1, 1), tag::flat<tag::dense>(), tag::flat<tag::multi_vector>(), tag::flat<tag::multi_vector>());

    // test(cmr(), tag::flat<tag::dense>()); // counterexample

    return 0; 
}

