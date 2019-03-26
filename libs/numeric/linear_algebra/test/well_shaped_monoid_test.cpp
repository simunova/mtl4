#include <iostream>
#include <cmath>


#ifdef __GXX_CONCEPTS__
#  include <concepts>
#else 
#  include <boost/numeric/linear_algebra/pseudo_concept.hpp>
#endif

struct vec
{
    vec (int s) : s(s) {}
    int s;
};

// Just to have another type
template <typename Value>
struct svec
{
    vec (int s) : s(s), p(new Value[s]) {}
    ~vec() { delete[](p); }
    int s;
  private:
    Value* p;
};







struct add {};

vec inline identity(add, vec v) { return v; }

concept WellShapedType<typename T>
{
    bool same_shape(T, T);

    axiom Reflexivity(T a) { same_shape(a, a) == true; }
    axiom Symmetry(T a, T b) 
    {
	if (same_shape(a, b))
	    same_shape(b, a) == true;
    }
    axiom Transitivity(T a, T b, T c) 
    {
	if (same_shape(a, b) && same_shape(b, c))
	    same_shape(a, c) == true;
    }
}


concept Monoid<typename Operation, typename Element>
{
    Element identity(Operation, Element);
}


concept WellShapedMonoid<typename Operation, typename Element>
  : Monoid<Operation, Element>
{
    requires WellShapedType<Element>;

    axiom ShapedIdentity(Operation op, Element x, Element y)
    {
	if (same_shape(x, y))
	    identity(op, x) == identity(op, y);
    }
}

concept_map WellShapedType<vec> {}
concept_map WellShapedMonoid<add, vec> {}


int main(int, char* [])  
{
    vec u(3), v(identity(add(), u));
   

    return 0;
}
