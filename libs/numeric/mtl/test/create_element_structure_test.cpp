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

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/operation/element_structure_algorithms.hpp>
#include <boost/numeric/itl/pc/matrix_algorithms.hpp>

int main()
{
#if 0    
    typedef mtl::compressed2D<double>    Matrix;
    typedef double value_type;
    
    Matrix A(mtl::io::matrix_market("/home/cornelius/projects/data/sysMat.mtx"));
    std::cout<< "A[0][0]=" << A[0][0] << "\n";
    std::cout<< "cols=" << num_cols(A) << "\n";
    std::cout<< "rows=" << num_rows(A) << "\n";
    
    mtl::mat::element_structure<value_type> es;
    std::string output("/home/cornelius/projects/data/sysMat_elem.mtx");
    imf::greedy_extract_element_structure(es, A, output);
    
     int size( es.get_total_vars() );
   
    mtl::dense_vector<value_type>              x(size, 1), b(size), ident(size); 
     iota(ident);
    mtl::compressed2D<double> B;
    assemble_compressed(es, B, ident);
    
    Matrix C(A-B);
    
    std::cout<< "norm(A-B)=" << frobenius_norm(C) << "\n";
#endif
    std::cout << "does nothing, more an example than a test\n";
    return 0;
}
