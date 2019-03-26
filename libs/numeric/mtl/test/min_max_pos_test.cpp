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

namespace std {

   template <typename T, typename U>
    inline ostream& operator<< (ostream& os, pair<T, U> const& p)
    {
        return os << '(' << p.first << ',' << p.second << ')';
    }
}
    
template <typename Matrix>
void init(Matrix& A)
{    
    A= 0.0;
    mtl::mat::inserter<Matrix>  ins(A);
    
    ins[0][1] << 3; ins[1][4] << 7;
    ins[2][3] << -2; ins[2][4] << 5;
}


template <typename Coll, typename Pos>
void test(Coll& coll, const char* name, Pos exp_min, Pos exp_max)
{
    cout << "\n" << name << " =\n" << coll; 

    Pos my_min= Pos(min_pos(coll));
    cout << "\nPosition of minimum is " << my_min << "\n";
    MTL_THROW_IF(my_min != exp_min, mtl::runtime_error("Minimum not at expected position"));

    Pos my_max= Pos(max_pos(coll));
    cout << "Position of maximum is " << my_max << "\n";
    MTL_THROW_IF(my_max != exp_max, mtl::runtime_error("Maximum not at expected position"));
}
 

int main(int, char**)
{
    using namespace mtl;

    dense_vector<double>  x(5);
    x= 1, 2, -3, 5, 4;
    dense2D<double>       A(3, 5);
    init(A);
    compressed2D<double>  B(3, 5);
    init(B);

    typedef mtl::Collection<dense2D<double> >::size_type size_type;
    std::pair<size_type, size_type> minp(2, 3), maxp(1, 4);

    test(x, "double vector", 2, 3);
    test(A, "dense matrix", minp, maxp);
    test(B, "sparse matrix", minp, maxp);

    return 0;
}
 














