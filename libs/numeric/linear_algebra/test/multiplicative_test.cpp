#include <iostream>
#include <cmath>
#include <concepts>

#include <boost/numeric/linear_algebra/operators.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>
#include <boost/numeric/linear_algebra/is_invertible.hpp>
#include <boost/numeric/linear_algebra/new_concepts.hpp>
#include <boost/numeric/linear_algebra/concept_maps.hpp>
#include <boost/numeric/linear_algebra/power.hpp>


template <std::HasMultiply T>
    requires std::Convertible<std::HasMultiply<T, T>::result_type, T>
          && std::Semiregular<T>
T f(const T& x, const T& y)
{
    std::cout << "\nMultiplicative magma.\n";
    return x * y;
}

template <math::MultiplicativeSemiGroup T>
    requires std::Convertible<std::HasMultiply<T, T>::result_type, T>
          && std::Semiregular<T>
T f(const T& x, const T& y)
{
    std::cout << "\nMultiplicative semi-group.\n";
    return x * y;
}

template <math::MultiplicativeMonoid T>
    requires std::Convertible<std::HasMultiply<T, T>::result_type, T>
          && std::Semiregular<T>
T f(const T& x, const T& y)
{
    std::cout << "\nMultiplicative monoid.\n";
    return x * y;
}

#if 0
template <math::MultiplicativePIMonoid T>
    requires std::Convertible<std::HasMultiply<T, T>::result_type, T>
          && std::Semiregular<T>
T f(const T& x, const T& y)
{
    std::cout << "\nMultiplicative partially invertible monoid.\n";
    return x * y;
}

template <math::MultiplicativeGroup T>
    requires std::Convertible<std::HasMultiply<T, T>::result_type, T>
          && std::Semiregular<T>
T f(const T& x, const T& y)
{
    std::cout << "\nMultiplicative group.\n";
    return x * y;
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

