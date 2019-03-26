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
#include <random>
#include <string>

#include <boost/numeric/mtl/mtl.hpp>


using namespace std;

int main(int argc, char** argv) 
{
    default_random_engine re(random_device{}());
    
    
    int n = 10;
    if (argc > 1)
        n = stoi(argv[1]);
        
    for (int s = 4; s < n; ++s) {
        cout << "s = " << s << endl;
        mtl::compressed2D<int> A(s, s), B(s, s);
        {
            mtl::mat::inserter<mtl::compressed2D<int>> ia(A, 5), ib(B, 0);
            uniform_int_distribution<> u(0, s-1);
            for (int i = 0; i < 4*s; ++i) {
                int r = u(re), c = u(re), v = u(re) + 1;
                ia[r][c] << v;
                ib[r][c] << v;
                // cout << "A[" << r << "][" << c << "] = " << v << endl;
            }
        }
        if (s < 8)  
            cout << "A =\n" << A << "B =\n" << B;
        mtl::compressed2D<int> D(A - B);
        if (one_norm(D) > 0) {
            cout << "Matrices are different!\n";
            return 1;
        }
    }

    return 0;
}
