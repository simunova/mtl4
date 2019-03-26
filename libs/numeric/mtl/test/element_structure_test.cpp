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
//
// Algorithm inspired by Nick Vannieuwenhoven, written by Cornelius Steinhardt

#define MTL_VPT_LEVEL 3

#include <iostream>
#include <boost/timer.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/mtl/io/read_el_matrix.hpp>



int main(int, char** argv)
{
    mtl::vampir_trace<9999> tracer;
    typedef double value_type;
       
    std::string program_dir= mtl::io::directory_name(argv[0]),
  	        matrix_file= mtl::io::join(program_dir, "../../mtl/test/matrix_market/square3.mtx");
// 		matrix_file= mtl::io::join(program_dir, "../../mtl/test/matrix_market/obstacle_small.mtx");
// 	        matrix_file= mtl::io::join(program_dir, "../../../../../branches/data/matrix_market/obstacle_q1q1_e64/obstacle_q1q1_e64_r00800.mtx");

    
    mtl::mat::element_structure<value_type> A;
    read_el_matrix(matrix_file, A);
    A.make_compact();

    const int size= A.get_total_vars();
    mtl::dense_vector<value_type>              x(size, 1), b(size), ident(size); 
    iota(ident);
    b= 7.0;
     
    for (int i= 0; i < 100; i++) {
	b+= A * x;
	if (size <= 100)
	    std::cout << "b = " << b << '\n';
    }

    if (size <= 100)
	std::cout << "b = " << b << '\n';

    return 0;
}
