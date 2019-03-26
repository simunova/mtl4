#include <iostream>
#include <cmath>


#ifdef __GXX_CONCEPTS__
#  include <concepts>
#else 
#  include <boost/numeric/linear_algebra/pseudo_concept.hpp>
#endif


namespace algebra {

    struct add {};

#if 0 // Doesn't work
    concept Monoid<typename Operation, typename Element>
    {
	Element identity<Operation, Element>();
    };

    template <typename Operation, typename Element>
    struct identity_functor {};

    template <typename Element>
    struct identity_functor<add, Element> 
    {
	Element operator()() { return 0; }
    };
    
    template <typename Operation, typename Element>
    Element identity()
    {
	return identity_functor<Operation, Element>()();
    }

    concept_map Monoid<add, int> {}
#endif

    concept Monoid<typename Operation, typename Element>
    {
	Element identity(Operation, Element);
    }

    concept StaticMonoid<typename Operation, typename Element>
    {
	Element static_identity();
	
	axiom Identity(Operation op, Element x)
	{
	    op(x, static_identity()) == x;
	    op(static_identity(), x) == x;
	}

	// Default implementation
	// Element identity(Operation, Element) { return static_identity(); }
    }
        
    template <typename Element>
    concept_map StaticMonoid<add, Element>
    {
	Element static_identity() { return Element(0); }
    }

#if 0
    template <typename Element>
    Element f_add(const Element& x, const Element& y) { return x + y; }

    template <typename Element>
    concept_map StaticMonoid<Element f_add<Element>(), Element>
    {
	Element static_identity() { return Element(0); }
    }
#endif
}



int main(int, char* [])  
{
    return 0;
}
