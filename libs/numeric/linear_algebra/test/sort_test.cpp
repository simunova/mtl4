#include <iostream>
#include <cmath>
#include <iterator>
#include <algorithm>


#ifdef __GXX_CONCEPTS__
#  include <concepts>
#else 
#  include <boost/numeric/linear_algebra/pseudo_concept.hpp>
#endif

concept PartialOrdering<typename Comparison, typename T>
{
    axiom Irreflexivity(Comparison cmp, T x)
    {
	!cmp(x, x);
    }

    axiom AntiSymmetry(Comparison cmp, T x, T y)
    {
	!(cmp(x, y) && cmp(y, x));
    }
	
    axiom Transitivity(Comparison cmp, T x, T y, T z)
    {
	if (cmp(x, y) && cmp(y, z))
	    cmp(x, z);
    }
}


concept StrictWeakOrdering<typename Comparison, typename T>
  : PartialOrdering<Comparison, T>
{
    axiom ComplementaryTransitivity(Comparison cmp, T x, T y, T z)
    {
	if (!cmp(x, y) && !cmp(y, z))
	    !cmp(x, z);
    }

    // Implies the following definition of equivalence classes
    bool equivalent(T x, T y)
    {
	return !cmp(x, y) && !cmp(y, x);
    }
}


concept TotalOrdering<typename Comparison, typename T>
  : StrictWeakOrdering<Comparison, T>
{
    axiom Trichotomy(Comparison cmp, T x, T y)
    {
	// Can we request that x == y is defined reasonably?
	cmp(x, y) || cmp(y, x) || x == y;
    }
}


// Same for operator

concept LessThanPartialOrdering<typename T>
{
    axiom Irreflexivity(T x)
    {
	!(x < x);
    }

    axiom AntiSymmetry(T x, T y)
    {
	!(x < y && y < x);
    }
	
    axiom Transitivity(T x, T y, T z)
    {
	if (x < y && y < z)
	    x < z;
    }
}


concept LessThanStrictWeakOrdering<typename T>
  : LessThanPartialOrdering<T>
{
    axiom ComplementaryTransitivity(T x, T y, T z)
    {
	if (!(x < y) && !(y < x))
	    !(x < z);
    }

    // Implies the following definition of equivalence classes
    bool equivalent(T x, T y)
    {
	return !(x < y) && !(y < x);
    }
}


concept LessThanTotalOrdering<typename T>
  : LessThanStrictWeakOrdering<T>
{
    axiom Trichotomy(T x, T y)
    {
	// Can we request that x == y is defined reasonably?
	x < y || y < x || x == y;
    }
}




int main(int, char* [])  
{
    int numbers[] = {7, 8, 3, 77, 11};
    int n= sizeof numbers / sizeof(char*);
    std::sort(numbers, numbers+n);

    std::copy(numbers, numbers+n, std::ostream_iterator<int>(std::cout, " ")); 
    std::cout << "\n";
    return 0;
}
