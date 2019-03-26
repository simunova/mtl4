#include <iostream>
#include <cmath>


#ifdef __GXX_CONCEPTS__
#  include <concepts>
#else 
#  include <boost/numeric/linear_algebra/pseudo_concept.hpp>
#endif


namespace algebra {

    concept Associative<typename Operation, typename Element>
      : std::Callable2<Operation, Element, Element>
    {
        axiom Associativity(Operation op, Element x, Element y, Element z)
        {
	    op(x, op(y, z)) == op(op(x, y), z); 
        }

        typename some_type;
        some_type some_function(Operation, Element);
    };


    auto concept SemiGroup<typename Operation, typename Element>
      : Associative<Operation, Element>
    {};

    concept Monoid<typename Operation, typename Element>
      : SemiGroup<Operation, Element> 
    {
    };

    auto concept Ring<typename AddOp, typename MultOp, typename Element>
      : Monoid<AddOp, Element>,
        SemiGroup<MultOp, Element>
    {};

    auto concept RingWithIdentity<typename AddOp, typename MultOp, typename Element>
      : Ring<AddOp, MultOp, Element>,
        Monoid<MultOp, Element>
    {};

}



int main(int, char* [])  
{
    return 0;
}
