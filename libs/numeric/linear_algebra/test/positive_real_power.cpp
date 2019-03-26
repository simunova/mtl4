#include <iostream>
#include <cmath>

#include <boost/numeric/linear_algebra/operators.hpp>
#include <boost/numeric/linear_algebra/concepts.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>
#include <boost/numeric/linear_algebra/is_invertible.hpp>

#include <libs/numeric/linear_algebra/test/algebraic_functions.hpp>
#include <libs/numeric/linear_algebra/test/power.hpp>
#include <libs/numeric/linear_algebra/test/positive_real.hpp>
#include <libs/numeric/linear_algebra/test/positive_real_power.hpp>

// User defined data types and operators

using mtl::positive_real;

// We assume that 0 and infinity is not to guarantee invertibility

# ifdef __GXX_CONCEPTS__
  namespace math { 
      concept_map SemiGroup< semigroup_mult, positive_real > {};
      concept_map Monoid< monoid_mult, positive_real > {};
      concept_map PartiallyInvertibleMonoid< pim_mult, positive_real > {};
      concept_map Group< group_mult, positive_real > {};
  }
# endif
 



int main(int, char* []) 
{
    using mtl::power;
    using math::mult;
 
    positive_real          value(1.1), zero(0.0);

    compute_power(value, 777, magma_mult(), "Magma");
    std::cout << '\n';

    compute_power(value, 777, semigroup_mult(), "SemiGroup");
    std::cout << '\n';

    compute_power(value, 777, monoid_mult(), "Monoid");
    compute_power(value, -777, monoid_mult(), "Monoid");
    std::cout << '\n';

    compute_power(value, 777, pim_mult(), "PIMonoid");
    compute_power(value, -777, pim_mult(), "PIMonoid");
    compute_power(zero, 777, pim_mult(), "PIMonoid");
    compute_power(zero, -777, pim_mult(), "PIMonoid");
    std::cout << '\n';

    compute_power(value, 777, group_mult(), "Group");
    compute_power(value, -777, group_mult(), "Group");
    compute_power(zero, 777, group_mult(), "Group");
    compute_power(zero, -777, group_mult(), "Group");
 
    return 0;
}

