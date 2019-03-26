// File: vector.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

template <typename Vector>
void fill_and_print(Vector& v, const char* name)
{
    // Set values in traditional way and print them
    v= 1.2, 3.4, 5.6;
    std::cout << name << " is " << v << "\n";
}

int main(int, char**)
{
#if defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_TEMPLATE_ALIAS)
    using namespace mtl;

    // Regular vector
    vector<double>                                  v1(3);
    fill_and_print(v1, "v1");

    // Row vector
    vector<double, row_major>                       v2(3);
    fill_and_print(v2, "v2");

    // Fixed-size vector
    vector<double, dim<3> >                         v3;
    fill_and_print(v3, "v3");
#endif
    return 0;
}

