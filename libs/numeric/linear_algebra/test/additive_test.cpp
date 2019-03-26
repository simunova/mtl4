#include <iostream>
#include <cmath>
#include <concepts>

#include <boost/numeric/linear_algebra/operators.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>
#include <boost/numeric/linear_algebra/is_invertible.hpp>
#include <boost/numeric/linear_algebra/new_concepts.hpp>
#include <boost/numeric/linear_algebra/concept_maps.hpp>
#include <boost/numeric/linear_algebra/power.hpp>


template <std::HasPlus T>
    requires std::Convertible<std::HasPlus<T, T>::result_type, T>
          && std::Semiregular<T>
T f(const T& x, const T& y)
{
    std::cout << "\nAdditive magma.\n";
    return x + y;
}

template <math::AdditiveSemiGroup T>
    requires std::Convertible<std::HasPlus<T, T>::result_type, T>
          && std::Semiregular<T>
T f(const T& x, const T& y)
{
    std::cout << "\nAdditive semi-group.\n";
    return x + y;
}

template <math::AdditiveMonoid T>
    requires std::Convertible<std::HasPlus<T, T>::result_type, T>
          && std::Semiregular<T>
T f(const T& x, const T& y)
{
    std::cout << "\nAdditive monoid.\n";
    return x + y;
}

#if 0
template <math::AdditivePIMonoid T>
    requires std::Convertible<std::HasPlus<T, T>::result_type, T>
          && std::Semiregular<T>
T f(const T& x, const T& y)
{
    std::cout << "\nAdditive partially invertible monoid.\n";
    return x + y;
}

template <math::AdditiveGroup T>
    requires std::Convertible<std::HasPlus<T, T>::result_type, T>
          && std::Semiregular<T>
T f(const T& x, const T& y)
{
    std::cout << "\nAdditive group.\n";
    return x + y;
}
#endif

int main(int, char* []) 
{
    std::cout << "f(3, 4) " << f(3, 4) << '\n';
    std::cout << "f(3l, 4l) " << f(3l, 4l) << '\n';
    std::cout << "f(3u, 4u) " << f(3u, 4u) << '\n';
    std::cout << "f(3.0f, 4.0f) " << f(3.0f, 4.0f) << '\n';
    std::cout << "f(3.0, 4.0) " << f(3.0, 4.0) << '\n';

    return 0;
}

