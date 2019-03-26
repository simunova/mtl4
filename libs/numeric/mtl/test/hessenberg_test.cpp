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
#include <utility>
#include <cmath>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;  
   
template <typename Matrix>
void init(Matrix& A)
{    
    A= 0.0;
    mtl::mat::inserter<Matrix>  ins(A);
    
    ins[0][1] << 3;  ins[1][4] << 7; ins[0][0] << 1; ins[4][4] << 17;
    ins[2][3] << -2; ins[2][4] << 5; ins[4][0] << 2; ins[4][1] <<  3;
    ins[3][2] << 4;
}


template <typename Coll>
void test(Coll& coll, const char* name)
{
    cout << "\n" << name << " =\n" << coll << "\n";
    
    Coll E(coll), F(coll); F=0;
    E= hessenberg(coll);
    F=tril(E,-2);
    std::cout<< "Hessenberg=\n" << E << "\n";
    std::cout<< "triu=\n" << F << "\n";
  
    MTL_THROW_IF(one_norm(F)>0.000004, mtl::runtime_error("No Hessenberg-Form"));
    F=extract_householder_hessenberg(coll);
    std::cout<< "extract_householder_hessenberg=\n" << F << "\n";
    F=extract_hessenberg(coll);
    std::cout<< "extract_hessenberg=\n" << F << "\n";
    F=householder_hessenberg(coll);
    std::cout<< "householder_hessenberg=\n" << F << "\n";
    F=hessenberg_factors(coll);
    std::cout<< "hessenberg_factors=\n" << F << "\n";
    //      F=hessenberg_q(coll);
    //      std::cout<< "hessenberg_q=\n" << F << "\n";
}
 

int main(int, char**)
{
    using namespace mtl;

    dense2D<double>       A(5, 5);
    init(A);
    compressed2D<double>  B(5, 5);
    init(B);

    test(A, "dense matrix");
//      test(B, "sparse matrix");  // TODO sparse Form

    return 0;
}
 














