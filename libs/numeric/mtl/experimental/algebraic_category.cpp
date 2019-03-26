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
#include <boost/test/minimal.hpp>

#include <boost/numeric/mtl/recursion/bit_masking.hpp>
#include <boost/numeric/mtl/recursion/predefined_masks.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp> 
#include <boost/numeric/mtl/matrix/compressed2D.hpp> 
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>


using namespace mtl;
using namespace std;  

template <typename T>
void dispatch(const T& x, tag::scalar, const char* name)
{
    cout << "dispatch " << name << " to scalar\n";
}

template <typename T>
void dispatch(const T& x, tag::vector, const char* name)
{
    cout << "dispatch " << name << " to vector\n";
}

template <typename T>
void dispatch(const T& x, tag::matrix, const char* name)
{ 
    cout << "dispatch " << name << " to matrix\n";
}


template <typename T>
void test(const T& x, const char* name)
{
    dispatch(x, typename traits::algebraic_category<T>::type(), name);
}

struct none {};


namespace mtl { namespace traits {

    template <typename Value>
    struct algebraic_category<Value[]>
    {
	typedef tag::vector type;
    };
}}

template <typename Value>
struct my_category
{
    typedef tag::scalar type;
};


template <typename Value, unsigned Size>
struct my_category<Value[Size]>
{
    typedef tag::vector type;
};

template <typename Value>
struct my_category<Value*>
{
    typedef tag::vector type;
};

template <typename Value, unsigned Rows, unsigned Cols>
struct my_category<Value[Rows][Cols]>
{
    typedef tag::matrix type;
};

template <typename> struct algebra_size;
template <typename> struct static_num_rows;
template <typename> struct static_num_cols;


template <typename T>
void test2(T const& x)
{
    cout << typeid(T).name() << "  " << typeid(typename my_category<T>::type).name() << "\n";
    algebra_size<T>()(x);
    cout << "static_num_rows = " << static_num_rows<T>::value << " static_num_cols = " << static_num_cols<T>::value << "\n\n";
}


template <typename T>
struct test3
{
    void operator()(T x)
    {
	cout << typeid(T).name() << "  " << typeid(typename my_category<T>::type).name() << "\n";
    }
};

template <typename T>
struct algebra_size 
{
    void operator()(const T& x)
    {
	cout << "Scalar\n";
    }
};


template <typename Value, unsigned Rows>
struct algebra_size<Value[Rows]>
{
    typedef Value  array_type[Rows];

    void operator()(const array_type& x)
    {
	cout << "size() = " << sizeof(x) / sizeof(x[0]) << "\n"; 
    }
};

template <typename Value, unsigned Rows, unsigned Cols>
struct algebra_size<Value[Rows][Cols]>
{
    typedef Value  array_type[Rows][Cols];

    void operator()(const array_type& x)
    {
	cout << "num_rows() = " << sizeof(x) / sizeof(x[0]) << " num_cols() = " << sizeof(x[0]) / sizeof(x[0][0]) << "\n"; 
    }
};

template <typename T>
struct static_num_rows 
{
    static const unsigned long value = 1;
};

template <typename T>
struct static_num_cols 
{
    static const unsigned long value = 1;
};



template <typename Value, unsigned Rows>
struct static_num_rows<Value[Rows]>
{
    static const unsigned long value = Rows;
};

template <typename Value, unsigned Rows>
struct static_num_cols<Value[Rows]>
{
    static const unsigned long value = 1;
};



template <typename Value, unsigned Rows, unsigned Cols>
struct static_num_rows<Value[Rows][Cols]>
{
    static const unsigned long value = Rows;
};

template <typename Value, unsigned Rows, unsigned Cols>
struct static_num_cols<Value[Rows][Cols]>
{
    static const unsigned long value = Cols;
};



int test_main(int argc, char* argv[])
{
    const int size= 5;
    double d;
    int    i;
    none   n;

    dense_vector<float>                                  u(size);
    dense_vector<float, mtl::parameters<row_major> >  ur(size);

    dense2D<double>                                      dr(size, size);
    dense2D<double, mat::parameters<col_major> >      dc(size, size);
    morton_dense<double, recursion::morton_z_mask>       mzd(size, size);
    morton_dense<double, recursion::doppled_2_row_mask>  d2r(size, size);
    compressed2D<double>                                 cr(size, size);
    compressed2D<double, mat::parameters<col_major> > cc(size, size);

    dense2D<complex<double> >                            drc(size, size);
    compressed2D<complex<double> >                       crc(size, size);

    mat::scaled_view<double, dense2D<double> >        scaled_matrix(2.0, dr);
    int   array[7];
    int   array2D[7][6];
    int   array3D[7][6][5];

    test2(i);
    test2(array);
    test2(array2D);
    test2(array3D);

    test3<int>()(i);
    test3<int[7]>()(array);
    test3<int[7][6]>()(array2D);
    test3<int[7][6][5]>()(array3D);

    algebra_size<int[7][6]>()(array2D);
    algebra_size<int[7][6][5]>()(array3D);

    return 0;

    test(d, "double");
    test(i, "int");
    test(n, "unknown type");

    test(u, "dense (column) vector");
    test(ur, "dense row vector");

    test(dr, "Dense row major");
    test(dc, "Dense column major");
    test(mzd, "Morton Z-order");
    test(d2r, "Hybrid 2 row-major");
    test(cr, "Compressed row major");
    test(cc, "Compressed column major");
    test(drc, "Dense row major complex");
    test(crc, "Compressed row major complex");
    test(scaled_matrix, "Scaled_view of dense row major");
    test(array, "array");
    

    return 0;
}
