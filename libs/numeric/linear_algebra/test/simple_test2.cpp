#include <iostream>
#include <cmath>


#  include <concepts>

namespace math {

// ==================================
// Classification of Arithmetic Types
// ==================================

// In addtion to std::Integral
concept Float<typename T> 
  : std::DefaultConstructible<T>, std::CopyConstructible<T>,
    std::EqualityComparable<T>
{
  T operator+(T);
  T operator+(T, T);
  T& operator+=(T&, T);
  T operator-(T, T);
  T operator-(T);
  T& operator-=(T&, T);
  T operator*(T, T);
  T& operator*=(T&, T);
  T operator/(T, T);
  T& operator/=(T&, T);


  //requires std::Assignable<T>, std::SameType<std::Assignable<T>::result_type, T&>;
}

concept_map Float<float> {}
concept_map Float<double> {}
concept_map Float<long double> {}


} // namespace math


int main(int, char* [])  
{

    return 0;
}
