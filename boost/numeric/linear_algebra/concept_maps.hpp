// $COPYRIGHT$

#ifndef MATH_CONCEPT_MAPS_INCLUDE
#define MATH_CONCEPT_MAPS_INCLUDE

#include <boost/numeric/linear_algebra/intrinsic_concept_maps.hpp>
#include <boost/numeric/linear_algebra/new_concepts.hpp>


namespace math {

#if 0 // Why this doesn't work???
    template <typename T>
        requires IntrinsicUnsignedIntegral<T>
    concept_map UnsignedIntegral<T> {}

    template <typename T>
        requires IntrinsicSignedIntegral<T>
    concept_map SignedIntegral<T> {}
#endif
    
    concept_map UnsignedIntegral<unsigned int> {}
    concept_map AdditiveMonoid<unsigned int> {}

    concept_map SignedIntegral<int> {}
    concept_map AdditiveSemiGroup<int> {}

    concept_map Commutative< add<int>, int > {}
    concept_map Monoid< add<int>, int > {}
    concept_map AbelianGroup< add<int>, int > {}

    concept_map Commutative< mult<int>, int > {}
    concept_map Monoid< mult<int>, int > {}
    concept_map PIMonoid< mult<int>, int > {}

    concept_map Commutative< min<int>, int > {}
    concept_map Monoid< min<int>, int > {}

    concept_map Commutative< max<int>, int > {}
    //concept_map Monoid< max<int>, int > {}


    concept_map Commutative< add<float>, float > {}
    concept_map Monoid< add<float>, float > {}
    concept_map AbelianGroup< add<float>, float > {}

    concept_map Commutative< mult<float>, float > {}
    concept_map Monoid< mult<float>, float > {}
    concept_map PIMonoid< mult<float>, float > {}

    concept_map AdditiveMonoid< float > {}
    // concept_map AdditiveAbelianGroup< float > {}

    // concept_map OperatorField<float> {}

    concept_map Commutative< min<float>, float > {}
    concept_map Monoid< min<float>, float > {}

    concept_map Commutative< max<float>, float > {}
    concept_map Monoid< max<float>, float > {}



    concept_map Commutative< add<double>, double > {}
    concept_map Monoid< add<double>, double > {}
    concept_map AbelianGroup< add<double>, double > {}

    concept_map Commutative< mult<double>, double > {}
    concept_map Monoid< mult<double>, double > {}
    concept_map PIMonoid< mult<double>, double > {}

    concept_map AdditiveMonoid< double > {}
    // concept_map AdditiveAbelianGroup< double > {}

    // concept_map OperatorField<double> {}

    concept_map Commutative< min<double>, double > {}
    concept_map Monoid< min<double>, double > {}

    concept_map Commutative< max<double>, double > {}
    concept_map Monoid< max<double>, double > {}



} // namespace math

#endif // MATH_CONCEPT_MAPS_INCLUDE
