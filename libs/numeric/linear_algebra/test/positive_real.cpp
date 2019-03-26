#include <iostream>
#include <cmath>

#include <boost/numeric/linear_algebra/operators.hpp>
#include <boost/numeric/linear_algebra/concepts.hpp>

#include <libs/numeric/linear_algebra/test/algebraic_functions.hpp>
#include <libs/numeric/linear_algebra/test/positive_real.hpp>

using mtl::positive_real;
 
namespace math {
    concept_map AdditiveMagma<mtl::positive_real> {};
}


int main(int, char* [])  
{
    positive_real                a0(0.0), a2(2.0), a3(3.0), a4(4.0), a5(5.0);
    math::add<positive_real>     my_add;
  
    std::cout << "equal_results(a2,a5,  a3,a4, my_add) " 
	      << mtl::equal_results(a2,a5,  a3,a4, my_add)  << std::endl;
    std::cout << "equal_results(a2,a4,  a3,a4, my_add) " 
	      << mtl::equal_results(a2,a4,  a3,a4, my_add)  << std::endl;

    std::cout << "equal_add_results(a2,a5,  a3,a4) " 
	      << mtl::equal_add_results(a2,a5,  a3,a4)  << std::endl;
    std::cout << "equal_add_results(a2,a4,  a3,a4) " 
	      << mtl::equal_add_results(a2,a4,  a3,a4)  << std::endl;

    return 0;
}

