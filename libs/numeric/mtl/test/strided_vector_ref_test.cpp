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
#include <boost/numeric/mtl/mtl.hpp>


using namespace std;  
    
template <typename Vector>
void one_d_iteration(char const* name, const Vector & vector, size_t check_index, typename Vector::value_type check)
{
    namespace traits = mtl::traits;
    typename traits::index<Vector>::type                               index(vector);
    typename traits::const_value<Vector>::type                         const_value(vector); 
    typedef  mtl::tag::nz                                              tag;
    typedef typename traits::range_generator<tag, Vector>::type        cursor_type;
    typedef typename traits::range_generator<tag, Vector>::complexity  complexity;

    cout << name << "\nElements: " << complexity() << '\n';
    for (cursor_type cursor = mtl::begin<tag>(vector), cend = mtl::end<tag>(vector); cursor != cend; ++cursor) {
	cout << "vector[" << index(*cursor) << "] = " 
	     << const_value(*cursor) << '\n';
	if (index(*cursor) == check_index 
	    && const_value(*cursor) != check) 
	    throw "wrong check value";
    }
}
    

template <typename Vector>
void test(Vector& v, const char* name)
{
    using mtl::sum; using mtl::product;

    typedef mtl::tag::iter::all iall;
    typedef typename mtl::traits::range_generator<iall, Vector>::type Iter;

    std::cout << "\n" << name << "  --- v = " << v << "\n"; std::cout.flush();

    for (Iter iter(mtl::begin<iall>(v)), iend(mtl::end<iall>(v)); iter != iend; ++iter)
	cout << *iter << ", "; 
    cout << "\n";

    one_d_iteration(name, v, 2, 8.0f);

    std::cout << "one_norm(v) = " << one_norm(v) << "\n"; std::cout.flush();
    MTL_THROW_IF(one_norm(v) != 15.0, mtl::runtime_error("one_norm wrong"));

    std::cout << "two_norm(v) = " << two_norm(v) << "\n"; std::cout.flush();
    MTL_THROW_IF(two_norm(v) < 9.6436 || two_norm(v) > 9.6437, mtl::runtime_error("two_norm wrong"));

    std::cout << "infinity_norm(v) = " << infinity_norm(v) << "\n"; std::cout.flush();
    MTL_THROW_IF(infinity_norm(v) != 8.0, mtl::runtime_error("infinity_norm wrong"));

    std::cout << "sum(v) = " << sum(v) << "\n"; std::cout.flush();
    MTL_THROW_IF(sum(v) != 15.0, mtl::runtime_error("sum wrong"));

    std::cout << "product(v) = " << product(v) << "\n"; std::cout.flush();
    MTL_THROW_IF(product(v) != 80.0, mtl::runtime_error("product wrong"));

    mtl::dense_vector<float>   u(3), w(3);
    u= 3.0; w= 5.0;

    std::cout << "u= v + w:\n"; std::cout.flush();
    u= v + w;
    cout << "u: " << u << "\n"; std::cout.flush();
    MTL_THROW_IF(u[0] != 7.0, mtl::runtime_error("wrong"));

    std::cout << "u= v + w + v + w:\n"; std::cout.flush();
    u= v + w + v + w;
    cout << "u: " << u << "\n"; std::cout.flush();
    MTL_THROW_IF(u[0] != 14.0, mtl::runtime_error("wrong"));

    std::cout << "u= v + (w= v + v);:\n"; std::cout.flush();
    u= v + (w= v + v);
    cout << "u: " << u << ", w: " << w << "\n"; std::cout.flush();
    MTL_THROW_IF(w[0] != 4.0, mtl::runtime_error("w wrong"));
    MTL_THROW_IF(u[0] != 6.0, mtl::runtime_error("u wrong"));

    std::cout << name << "  --- u+= dot<12>(v, w) * w;:\n"; std::cout.flush();
    u+= mtl::dot<12>(v, w) * w;
    cout << "u: " << u << ", v: " << v << ", w: " << w << "\n"; std::cout.flush();
    MTL_THROW_IF(u[0] != 750.0, mtl::runtime_error("u wrong"));

    Vector q(sub_vector(v, 1, 4));
    MTL_THROW_IF(q[1] != 8.f, mtl::runtime_error("Wrong value in q"));
    MTL_THROW_IF(size(q) != 2, mtl::runtime_error("Wrong size of q"));
    std::cout << "sub_vector(v, 1, 4) == " << q << "\n";
    

    using mtl::irange;
    Vector r(v[irange(1, 4)]);
    MTL_THROW_IF(r[1] != 8.f, mtl::runtime_error("Wrong value in r"));
    MTL_THROW_IF(size(r) != 2, mtl::runtime_error("Wrong size of r"));
    std::cout << "v[irange(1, 4)] == " << r << "\n";
    
}
 
template <typename Vector>
void test2(Vector& v)
{
#if 0
    typedef typename mtl::traits::category<Vector>::type cat;
    typedef typename mtl::traits::is_vector<Vector>::type isv;
    typedef typename mtl::Collection<Vector>::value_type value_type;
    typedef typename mtl::traits::category<value_type>::type vcat;
    typedef typename mtl::ashape::ashape<value_type>::type ashape;

    std::cout << "Category of v is " << typeid(cat).name() << '\n';
    std::cout << "is vector " << typeid(isv).name() << '\n';
    std::cout << "value type is " << typeid(value_type).name() << '\n';
    std::cout << "Category of value is " << typeid(vcat).name() << '\n';
    std::cout << "ashape of value is " << typeid(ashape).name() << '\n';
#endif
    set_to_zero(v);
}

int main(int, char**)
{
    mtl::dense_vector<float> v(3);
    v= 2, 5, 8;
    test(v, "Reference");

    mtl::dense2D<float> A(3, 3);
    A= 1, 2, 3,
       4, 5, 6, 
       7, 8, 9;
    
    mtl::vec::strided_vector_ref<float> stref(A[mtl::iall][2]);
    test2(v); 
    test2(stref); 

    mtl::vec::strided_vector_ref<float> x(3, &A[0][1], 3);

    const mtl::dense2D<float> B(A);
    mtl::vec::strided_vector_ref<const float> xc(3, &B[0][1], 3);

    test(x, "test float");
    test(xc, "test const float");

    return 0;
}
 












