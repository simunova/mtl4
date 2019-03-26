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
#include <boost/timer.hpp>


using namespace std;
using namespace mtl;




typedef compressed2D<double> sp_mat;


void assemble(sp_mat& A, double val)
{
  mat::inserter<sp_mat, update_plus<double> > ins(A, 3);

  double array[][3]= {{-2*val, val, val},
		      {val, -2*val, val},
		      {val, val, -2*val}};

  dense2D<double> block(array);
  dense_vector<int> cols(3);
  dense_vector<int> rows(3);

  int N= num_rows(A); // A is N-by-N
  for(int k=0; k<N-2; k+=3)    // nice case
  //for(int k=N-3; k>=0; k-=3) // out-of-order case
    {
      rows[0] = k; rows[1] = k+1; rows[2] = k+2;
      cols[0] = 0; cols[1] = k;   cols[2] = N-1;

      ins << element_matrix(block, rows, cols);
    }
}



int main(int argc, char* argv[])
{
    sp_mat B(9, 9);
    assemble(B, 1.0);
    cout << "Small assembled matrix\n" << B << "\n";

    const int size= 900000; 
    sp_mat A(size, size);

    boost::timer atime;
    assemble(A, 1.0);	
    cout << "Assemble time = " << atime.elapsed() << ", " 
	 << 3*size / atime.elapsed() << " elements/s\n";

    A*= 0.0;
    boost::timer rtime;
    assemble(A, 1.0);	
    cout << "Reassemble time = " << rtime.elapsed() << ", " 
	 << 3*size / rtime.elapsed() << " elements/s\n";

    return 0;
}
