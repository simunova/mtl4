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

#  define MTL_VPT_LEVEL 5

#include <iostream>
#include <limits>

#include <boost/timer.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <boost/numeric/mtl/io/read_el_matrix.hpp>
#include <boost/numeric/itl/pc/matrix_algorithms.hpp>
#include <boost/numeric/itl/pc/imf_preconditioner.hpp>
#include <boost/numeric/itl/pc/imf_algorithms.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#ifdef MTL_HAS_VPT
  #include <vt_user.h> 
#endif 

template< class ElementStructure >
void setup(ElementStructure& A, int lofi)
{
    typedef double value_type; 
      
    int size( A.get_total_vars() );
   
    mtl::dense_vector<value_type>              x(size, 1), b(size), ident(size); 
     iota(ident);

     std::cout<< "read ok, start with precond\n";
    boost::timer factorization;
    itl::pc::imf_preconditioner<value_type> precond(A, lofi);
    double ftime= factorization.elapsed();
    std::cout<< "imf ready\n";
    
//     for (int i= 0; i < size; i++)
//       x[i]= std::numeric_limits<value_type>::quiet_NaN();

    std::cout<< "size(rhs2)=" << num_rows(b) << "\n";
    mtl::compressed2D<double> B;
    assemble_compressed(A, B, ident);
    std::cout << "NNZ == " << B.nnz() << "\n";
    
    b= B * x;

#if 0
	mtl::io::tout << "------------------------------- STATISTICS -------------------------------" << std::endl;
	int rows = num_rows(*master_mat);
	int cols = num_cols(*master_mat);
 	int nnz = (*master_mat).nnz();
	mtl::io::tout << "Dimensions: " << rows << " x " << cols << std::endl;
	mtl::io::tout << "Non-zeros: " << nnz << std::endl;
	mtl::io::tout << "Sparsity (%): " << ((double(nnz) / rows) / cols) << std::endl;
	mtl::io::tout << "Avg nnz/row: " << double(nnz) / rows << std::endl;
	mtl::io::tout << std::endl;
	mtl::io::tout << "Elements: " << es.get_total_elements() << std::endl;
	mtl::io::tout << "Variables: " << es.get_total_vars() << std::endl;
	mtl::io::tout << "--------------------------------------------------------------------------" << std::endl;
// calculate eigenvalues
	mtl::dense2D<value_type> E(size,size),A(*master_mat);
	for(int i=0; i<size;i++){
	  mtl::dense_vector<value_type> tmp(A[mtl::irange(0, mtl::imax)][i]);
	  E[mtl::irange(0, mtl::imax)][i] = precond.solve(tmp);
	}
	mtl::io::tout<< "E=\n"<<E <<"\n";
#endif

	std::cout<< "start solver\n";
    itl::cyclic_iteration<value_type>          iter(b, size, 1.e-8, 0.0, 5);
    x= 0;
    boost::timer solver;
     bicgstab(B, x, b, precond, iter);
     std::cout << "Factorization took " << ftime << "s, solution took " << solver.elapsed() << "s\n";
}

int main(int, char** argv)
{
    mtl::vampir_trace<9999> tracer;
    typedef double value_type;
       
    std::string program_dir= mtl::io::directory_name(argv[0]),
  	        matrix_file= mtl::io::join(program_dir, "../../mtl/test/matrix_market/square3.mtx");
//  	        matrix_file= mtl::io::join(program_dir, "../../../../../data/sysMat_elem.mtx");
// 		matrix_file= mtl::io::join(program_dir, "../../mtl/test/matrix_market/obstacle_small.mtx");
// 	        matrix_file= mtl::io::join(program_dir, "../../../../../data/matrix_market/obstacle_q1q1_e64/obstacle_q1q1_e64_r00800.mtx");

    mtl::mat::element_structure<value_type> A;
    read_el_matrix(matrix_file, A);
	
    setup(A, 5);
    return 0;
}
