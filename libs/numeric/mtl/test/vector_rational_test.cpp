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



#include <boost/numeric/mtl/mtl.hpp>
# include "boost/rational.hpp"
# include "boost/range.hpp"



typedef boost::rational<long>                              t_Q;
typedef mtl::dense_vector<t_Q, mtl::vec::parameters<> > t_dVecQ;

//! Struct containing data about the refinment of a face.
struct st_test
{
    t_dVecQ vecQ;

    //! Constructor for refinement of a triangle.
    st_test( const t_dVecQ& vecQ_new )
    {
        vecQ = vecQ_new;
        std::cout << "vecQ: " << vecQ << "\n";       // uses size -> ambiguity with size in boost/range
        std::cout << mtl::size(vecQ_new) << "\n";       
    }
};

// using test as function name causes ambiguities  on some compilers
// (type_traits/has_new_operator.hpp:24: error: ‘template<class U, U x> struct boost::detail::test' is not a function)
void assign_test( const t_dVecQ& vecQ ) 
{
    t_dVecQ vecQ_temp;
    vecQ_temp = vecQ;
    std::cout << "vecQ_temp: " << vecQ_temp << "\n";
}

void test2()
{
   t_dVecQ vQ0(2,3);
    t_dVecQ vQ1(2,2);
    t_dVecQ vQ2(2,1);
    std::cout << "vQ0: " << vQ0 << "\n";
    std::cout << "vQ1: " << vQ1 << "\n";
    std::cout << "size(vQ1): " << mtl::size(vQ1) << "\n";
    vQ0 = vQ1 + vQ2;
    vQ0 = vQ1 - vQ2;
    vQ0 += vQ1;
    vQ0 -= vQ1;
    std::cout << "vQ0: " << vQ0 << "\n";
    t_dVecQ vQ3( vQ0 - vQ1 );
    std::cout << "vQ2: " << vQ2 << "\n";
    std::cout << "size(vQ2): " << mtl::size(vQ2) << "\n";
}


int main(int , char**)
{
    t_Q Q0(1,3);
    std::cout << "Q0: " << Q0 << "\n";

    t_dVecQ vecQ0(2,Q0);
    std::cout << "vecQ0: " << vecQ0 << "\n";

    t_Q Q1(1,2);
    std::cout << "Q1: " << Q1 << "\n";
    t_dVecQ vecQ1(3,Q1);
    std::cout << "vecQ1: " << vecQ1 << "\n";   
    vecQ1[mtl::irange(2)] = vecQ0;
    std::cout << "vecQ1: " << vecQ1 << "\n";   

    assign_test(vecQ1);
    st_test test(vecQ1);

    // std::cout << "size(vecQ1): " << mtl::size(vecQ1) << "\n";

    test2();

    return 0;
}
