#include <iostream>
#include <cmath>

#include <boost/numeric/linear_algebra/operators.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>
#include <boost/numeric/linear_algebra/is_invertible.hpp>
#include <boost/numeric/linear_algebra/new_concepts.hpp>
//#include <boost/numeric/linear_algebra/concept_maps.hpp>
#include <boost/numeric/linear_algebra/power.hpp>

#if 0 // only for testing other concepts
namespace math {
    concept_map AbelianGroup< add<float>, float > {}
    concept_map Monoid< mult<float>, float > {}
    concept_map SignedIntegral<int> {}
} // namespace math
#endif

int main(int, char* []) 
{
    using math::power; using math::mult; using math::add;
    float a= 3.14;

    std::cout << "power(a, 5, mult<float>) " << power(a, 5, mult<float>()) << '\n';
#if 0
    std::cout << "power(a, 5, add<float>) " << power(a, 5, add<float>()) << '\n';
    std::cout << "power(a, 0, add<float>) " << power(a, 0, add<float>()) << '\n';
#endif

    return 0;
}

