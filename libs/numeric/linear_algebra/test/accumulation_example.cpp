
#include <libs/numeric/linear_algebra/test/accumulation.hpp>

#include <boost/timer.hpp>


concept AccurateArithmetic<typename T> {}

template <typename T>
    requires std::Integral<T>
concept_map AccurateArithmetic<T> {}


concept TolerateRoundingErrors<typename T> {}

#define CONSIDER_FLOAT_ROUNDING_ERRORS

# ifndef CONSIDER_FLOAT_ROUNDING_ERRORS
    template <typename T>
        requires math::Float<T>
    concept_map TolerateRoundingErrors<T> {}
# endif

concept_map TolerateRoundingErrors<double> {}

concept SelectiveOperation<typename Operation, typename Element> {}
 // : std::CopyConstructible<std::vector<std::vector<Element> > >


template <typename Element>
concept_map SelectiveOperation<math::min<Element>, Element> {}

template <typename Element>
concept_map SelectiveOperation<math::max<Element>, Element> {}

concept RegularReduction<typename Operation, typename Element> {} // : std::Integral<float> {}

template <typename Operation, typename Element>
  requires AccurateArithmetic<Element>
concept_map RegularReduction<Operation, Element> {}

template <typename Operation, typename Element>
  requires TolerateRoundingErrors<Element>
concept_map RegularReduction<Operation, Element> {}

template <typename Operation, typename Element>
  requires SelectiveOperation<Operation, Element> 
        && !AccurateArithmetic<Element> 
        && !TolerateRoundingErrors<Element>
concept_map RegularReduction<Operation, Element> {}




// ######################################################################

namespace mtl {


// Dispatching between simple and unrolled version
template <typename Iter, typename Value, typename Op>
  requires std::ForwardIterator<Iter> 
                  && std::Convertible<Value, std::ForwardIterator<Iter>::value_type>
                  && math::Magma<Op, std::ForwardIterator<Iter>::value_type>
                  && RegularReduction<Op, std::ForwardIterator<Iter>::value_type>
typename std::ForwardIterator<Iter>::value_type 
inline my_accumulate(Iter first, Iter last, Value init, Op op)
{
    std::cout << "Simple accumulate\n";
    return mtl::accumulate_simple(first, last, init, op);
}


template <typename Iter, typename Value, typename Op>
    requires  std::RandomAccessIterator<Iter> 
	          && std::Convertible<Value, std::RandomAccessIterator<Iter>::value_type>
		  && math::Monoid<Op, std::RandomAccessIterator<Iter>::value_type> 
		  && math::Commutative<Op, std::RandomAccessIterator<Iter>::value_type> 
                  && RegularReduction<Op, std::RandomAccessIterator<Iter>::value_type>
typename std::RandomAccessIterator<Iter>::value_type 
inline my_accumulate(Iter first, Iter last, Value init, Op op)
{
    std::cout << "Unrolled accumulate\n";
    return mtl::accumulate_unrolled(first, last, init, op);
}

// Special Treatment
template <typename Iter, typename Value, typename Op>
  requires std::ForwardIterator<Iter> 
                  && std::Convertible<Value, std::ForwardIterator<Iter>::value_type>
                  && math::Magma<Op, std::ForwardIterator<Iter>::value_type>
typename std::ForwardIterator<Iter>::value_type 
inline my_accumulate(Iter first, Iter last, Value init, Op op)
{
    std::cout << "Special accumulate\n";
    return mtl::accumulate_simple(first, last, init, op);
}

} // namespace mtl

// ######################################################################

using math::identity; using math::add;

const int   array_size= 10;

template <typename Element>
void test_accumulate(const char* name)
{
    Element     array[array_size];
    for (int i= 0; i < array_size; i++) 
	array[i]= Element(i);

    std::cout << '\n' << name << '\n' << " Add: ";
    mtl::my_accumulate(&array[0], array+array_size, Element(0), math::add<Element>());
    std::cout << "Mult: ";
    mtl::my_accumulate(array, array+array_size, Element(1), math::mult<Element>());
    std::cout << " Min: ";
    mtl::my_accumulate(array, array+array_size, Element(1), math::min<Element>());
    std::cout << " Max: ";
    mtl::my_accumulate(array, array+array_size, Element(1), math::max<Element>());
}


int main(int, char* [])
{
    test_accumulate<int>("int");
    test_accumulate<float>("float");
    test_accumulate<double>("double");

    return 0;
}
