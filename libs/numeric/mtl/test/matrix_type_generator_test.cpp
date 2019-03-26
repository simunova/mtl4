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
#include <typeinfo>

#include <boost/type_traits/is_same.hpp>
#include <boost/numeric/mtl/mtl.hpp>

template <typename Matrix>
void compressed2D_check(const Matrix&)
{
    throw "Matrix should be compressed2D!\n";
}

template <typename Value, typename Parameter>
void compressed2D_check(const mtl::compressed2D<Value, Parameter>&) {}

template <typename Matrix>
void dense2D_check(const Matrix&)
{
    throw "Matrix should be dense2D!\n";
}

template <typename Value, typename Parameter>
void dense2D_check(const mtl::dense2D<Value, Parameter>&) {}



int main(int, char**)
{
    using namespace std;
    using namespace mtl;
    using mtl::io::tout;

#if defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_STATICASSERT) && defined(MTL_WITH_TEMPLATE_ALIAS)
    matrix<double, sparse> B;
    tout << "Type of B is " << typeid(B).name() << '\n';
    compressed2D_check(B);
    static_assert(traits::is_row_major<decltype(B)>::value, "Matrix must be row-major!");

    matrix<double, sparse, as_size_type<unsigned> > C;
    tout << "Type of C is " << typeid(C).name() << '\n';
    compressed2D_check(C);

    matrix<double, sparse, column_major> D;
    tout << "Type of D is " << typeid(D).name() << '\n';
    compressed2D_check(D);
    static_assert(!traits::is_static<decltype(D)>::value, "Matrix must not have static size!");
    static_assert(!traits::is_row_major<decltype(D)>::value, "Matrix must be column-major!");

    matrix<double, sparse, column_major, as_size_type<int>, dim<3, 5> > E;
    tout << "Type of E is " << typeid(E).name() << '\n';
    compressed2D_check(E);
    static_assert(traits::is_static<decltype(E)>::value, "Matrix must have static size!");
    static_assert(boost::is_same<mtl::Collection<decltype(E)>::size_type, int>::value, 
		  "Size type must be int!");



    matrix<double> F;
    tout << "Type of F is " << typeid(F).name() << '\n';
    dense2D_check(F);
    static_assert(traits::is_row_major<decltype(F)>::value, "Matrix must be row-major!");

    matrix<double, as_size_type<unsigned> > G;
    tout << "Type of G is " << typeid(G).name() << '\n';
    dense2D_check(G);

    matrix<double, column_major> I;
    tout << "Type of I is " << typeid(I).name() << '\n';
    dense2D_check(I);
    static_assert(!traits::is_static<decltype(I)>::value, "Matrix must not have static size!");
    static_assert(!traits::is_row_major<decltype(I)>::value, "Matrix must be column-major!");

    matrix<double, column_major, as_size_type<int>, dim<3, 5> > J;
    tout << "Type of J is " << typeid(J).name() << '\n';
    dense2D_check(J);
    static_assert(traits::is_static<decltype(J)>::value, "Matrix must have static size!");

    matrix<double, dim<3, 3> > K;
    tout << "Type of K is " << typeid(K).name() << '\n';
    dense2D_check(K);
    static_assert(traits::is_static<decltype(K)>::value, "Matrix must have static size!");

    matrix<double, banded>     L;
    tout << "Type of L is " << typeid(L).name() << '\n';
    static_assert(traits::is_sparse<decltype(L)>::value, "Matrix must be sparse!"); // Will be changed in the future
   
    matrix<double, coordinate>     M(7, 8);
    tout << "Type of M is " << typeid(M).name() << '\n';
    static_assert(traits::is_sparse<decltype(M)>::value, "Matrix must be sparse!"); 

    matrix<double, ellpack>     N;
    tout << "Type of N is " << typeid(N).name() << '\n';
    static_assert(traits::is_sparse<decltype(N)>::value, "Matrix must be sparse!"); 
   
    matrix<double, morton>     O;
    tout << "Type of O is " << typeid(O).name() << '\n';
    static_assert(!traits::is_sparse<decltype(O)>::value, "Matrix must be dense!"); 
    
    matrix<double, morton, dim<3, 5> >     P;
    tout << "Type of P is " << typeid(P).name() << '\n';
    static_assert(!traits::is_sparse<decltype(P)>::value, "Matrix must be dense!"); 
    


#else
    std::cout << "Test deactivated due to missing C++11 features.\n";
#endif

    return 0;
}
