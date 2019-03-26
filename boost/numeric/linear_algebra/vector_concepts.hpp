// Copyright 2006. Peter Gottschling, Matthias Troyer, Rolf Bonderer
// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#ifndef LA_VECTOR_CONCEPTS_INCLUDE
#define LA_VECTOR_CONCEPTS_INCLUDE


#include <boost/numeric/linear_algebra/concepts.hpp>
#include <boost/numeric/linear_algebra/ets_concepts.hpp>

#ifdef __GXX_CONCEPTS__
#  include <concepts>
#else 
#  include <boost/numeric/linear_algebra/pseudo_concept.hpp>
#endif


namespace math {  
  
/** @addtogroup Concepts
 *  @{
 */

#ifdef __GXX_CONCEPTS__
concept VectorSpace<typename Vector, typename Scalar = typename Vector::value_type>
: AdditiveAbelianGroup<Vector>
{
    requires Field<Scalar>;
    requires Multiplicable<Scalar, Vector>;
    requires MultiplicableWithAssign<Vector, Scalar>;
    requires DivisibleWithAssign<Vector, Scalar>;
  
    requires std::Assignable<Vector, Multiplicable<Scalar, Vector>::result_type>;
    requires std::Assignable<Vector, Multiplicable<Vector, Scalar>::result_type>;
    requires std::Assignable<Vector, Divisible<Vector, Scalar>::result_type>;
    
    // Associated types of Field<Scalar> and AdditiveAbelianGroup<Vector> collide
    // typename result_type = AdditiveAbelianGroup<Vector>::result_type;
    // typename assign_result_type = AdditiveAbelianGroup<Vector>::assign_result_type;

    axiom Distributivity(Vector v, Vector w, Scalar a, Scalar b)
    {
	a * (v + w) == a * v + a * w;
	(a + b) * v == a * v + b * v;
	// The following properties are implied by the above, Field and Abelian group
	// Can we be sure that compilers can deduce/interfere it?
	(v + w) * a == v * a + w * a;
	v * (a + b) == v * a + v * b;
    }
}
#else
    //! Concept VectorSpace
    /*!
	\param Vector   The the type of a vector or a collection 
        \param Scalar   The scalar over which the vector field is defined
        
	\par Requires:
	- Field < Scalar >;
	- Multiplicable <Scalar, Vector>;
	- MultiplicableWithAssign <Vector, Scalar>;
	- DivisibleWithAssign <Vector, Scalar>;
	- std::Assignable <Vector, Multiplicable<Scalar, Vector>::result_type>;
	- std::Assignable <Vector, Multiplicable<Vector, Scalar>::result_type>;
	- std::Assignable <Vector, Divisible<Vector, Scalar>::result_type>;

    */
    template <typename Vector, typename Scalar = typename Vector::value_type>
    struct VectorSpace
      : AdditiveAbelianGroup<Vector>
    {
	/// Invariant: Distributivity of scalars and vectors from left and from right
	axiom Distributivity(Vector v, Vector w, Scalar a, Scalar b)
	{
	    /// a * (v + w) == a * v + a * w;       // Scalar from left

	    /// Vector from right: (a + b) * v == a * v + b * v; 

	    /// Scalar from right: (v + w) * a == v * a + w * a; 

	    /// Vector from left:  v * (a + b) == v * a + v * b;
	}
    };
#endif

#ifdef __GXX_CONCEPTS__
concept Norm<typename N, typename Vector, 
	     typename Scalar = typename Vector::value_type>
  : std::Callable1<N, Vector>
{
    requires VectorSpace<Vector, Scalar>;
    requires RealMagnitude<Scalar>;
    typename magnitude_type = MagnitudeType<Scalar>::type;
    requires std::Convertible<magnitude_type, Scalar>;

    typename result_type_norm = std::Callable1<N, Vector>::result_type;
    requires std::Convertible<result_type_norm, RealMagnitude<Scalar>::magnitude_type>;
    requires std::Convertible<result_type_norm, Scalar>;

    // Version with function instead functor, as used by Rolf and Matthias
    // Axioms there defined without norm functor and concept has only 2 types
#if 0       
    typename result_type_norm; 
    result_type_norm norm(const Vector&);
    requires std::Convertible<result_type_norm, magnitude_type>;
    requires std::Convertible<result_type_norm, Scalar>;
#endif

    axiom Positivity(N norm, Vector v, magnitude_type ref)
    {
	norm(v) >= zero(ref);
    }

    // The following is covered by RealMagnitude
    // requires AbsApplicable<Scalar>;
    // requires std::Convertible<AbsApplicable<Scalar>::result_type, magnitude_type>;
    // requires Multiplicable<magnitude_type>;

    axiom PositiveHomogeneity(N norm, Vector v, Scalar a)
    {
	norm(a * v) == abs(a) * norm(v);
    }

    axiom TriangleInequality(N norm, Vector u, Vector v)
    {
	norm(u + v) <= norm(u) + norm(v);
    }
}
#else
    //! Concept Norm
    /*!
        Semantic requirements of a norm

	\param N        Norm functor
	\param Vector   The the type of a vector or a collection 
        \param Scalar   The scalar over which the vector field is defined
        
	\par Refinement of:
	- std::Callable1 <N, Vector>

	\par Associated types:
	- magnitude_type
	- result_type_norm

	\par Requires:
	- VectorSpace <Vector, Scalar>;
	- RealMagnitude < Scalar >;
	- std::Convertible <magnitude_type, Scalar>;
	- std::Convertible <result_type_norm, RealMagnitude<Scalar>::magnitude_type>;
	- std::Convertible <result_type_norm, Scalar>;

    */
template <typename N, typename Vector, 
	  typename Scalar = typename Vector::value_type>
struct Norm
  : std::Callable1<N, Vector>
{
    /// Associated type to represent real values in teh Field of scalar (with default)
    /** By default MagnitudeType<Scalar>::type */
    typedef associated_type magnitude_type;

    /// Associated type for result of norm functor
    /** Automatically detected */
    typedef associated_type result_type_norm;

    /// Invariant: norm of vector is larger than zero 
    axiom Positivity(N norm, Vector v, magnitude_type ref)
    {
	/// norm(v) >= zero(ref);
    }

    /// Invariant: positive homogeneity with scalar
    axiom PositiveHomogeneity(N norm, Vector v, Scalar a)
    {
	/// norm(a * v) == abs(a) * norm(v);
    }

    /// Invariant: triangle inequality
    axiom TriangleInequality(N norm, Vector u, Vector v)
    {
	/// norm(u + v) <= norm(u) + norm(v);
    }
};
#endif


#ifdef __GXX_CONCEPTS__
concept SemiNorm<typename N, typename Vector, 
		 typename Scalar = typename Vector::value_type>
  : Norm<N, Vector, Scalar>
{
    axiom PositiveDefiniteness(N norm, Vector v, magnitude_type ref)
    {
	if (norm(v) == zero(ref))
	    v == zero(v);
	if (v == zero(v))
	    norm(v) == zero(ref);
    }
}
#else
    //! Concept SemiNorm
    /*!
        Semantic requirements of a semi-norm

	\param N        Norm functor
	\param Vector   The the type of a vector or a collection 
        \param Scalar   The scalar over which the vector field is defined
        
	\par Refinement of:
	- Norm <N, Vector, Scalar>
    */
template <typename N, typename Vector, 
	  typename Scalar = typename Vector::value_type>
struct SemiNorm
  : Norm<N, Vector, Scalar>
{
    /// The norm of a vector is zero if and only if the vector is the zero vector
    axiom PositiveDefiniteness(N norm, Vector v, magnitude_type ref)
    {
	/// if (norm(v) == zero(ref)) v == zero(v);

	/// if (v == zero(v)) norm(v) == zero(ref);
    }
};
#endif

#ifdef __GXX_CONCEPTS__
concept BanachSpace<typename N, typename Vector, 
		    typename Scalar = typename Vector::value_type>
  : Norm<N, Vector, Scalar>,
    VectorSpace<Vector, Scalar>
{};
#else
    //! Concept BanachSpace
    /*!
        A Banach space is a vector space with a norm

	\param N        Norm functor
	\param Vector   The the type of a vector or a collection 
        \param Scalar   The scalar over which the vector field is defined
        
	\par Refinement of:
	- Norm <N, Vector, Scalar>
	- VectorSpace <Vector, Scalar>

	\note
	- The (expressible) requirements of Banach Space are already given in Norm.
	- The difference between the requirements is the completeness of the 
	  Banach space, i.e. that every Cauchy sequence w.r.t. norm(v-w) has a limit
	  in the space. Unfortunately, completeness is never satisfied for
	  finite precision arithmetic types.
	- Another subtle difference is that Norm is not a refinement of Vectorspace
    */
template <typename N, typename Vector, 
	  typename Scalar = typename Vector::value_type>
struct BanachSpace
  : Norm<N, Vector, Scalar>,
    VectorSpace<Vector, Scalar>
{};
#endif


#ifdef __GXX_CONCEPTS__
concept InnerProduct<typename I, typename Vector, 
		     typename Scalar = typename Vector::value_type>
  : std::Callable2<I, Vector, Vector>
{
    // Result of the inner product must be convertible to Scalar
    requires std::Convertible<std::Callable2<I, Vector, Vector>::result_type, Scalar>;

    // Let's try without this
    // requires ets::InnerProduct<I, Vector, Scalar>;

    requires HasConjugate<Scalar>;

    axiom ConjugateSymmetry(I inner, Vector v, Vector w)
    {
	inner(v, w) == conj(inner(w, v));
    }

    axiom SequiLinearity(I inner, Scalar a, Scalar b, Vector u, Vector v, Vector w)
    {
	inner(v, b * w) == b * inner(v, w);
	inner(u, v + w) == inner(u, v) + inner(u, w);
	// This implies the following (will compilers infere/deduce?)
	inner(a * v, w) == conj(a) * inner(v, w);
	inner(u + v, w) == inner(u, w) + inner(v, w);
    }

    requires RealMagnitude<Scalar>;
    typename magnitude_type = RealMagnitude<Scalar>::type;
    // requires FullLessThanComparable<magnitude_type>;

    axiom NonNegativity(I inner, Vector v, MagnitudeType<Scalar>::type magnitude)
    {
	// inner(v, v) == conj(inner(v, v)) implies inner(v, v) is real
	// ergo representable as magnitude type
	magnitude_type(inner(v, v)) >= zero(magnitude)
    }

    axiom NonDegeneracy(I inner, Vector v, Vector w, Scalar s)
    {
	if (v == zero(v))
	    inner(v, w) == zero(s);
	if (inner(v, w) == zero(s))
	    v == zero(v);
    }
};
#else
    //! Concept InnerProduct
    /*!
        Semantic requirements of a inner product

	\param I        The inner product functor
	\param Vector   The the type of a vector or a collection 
        \param Scalar   The scalar over which the vector field is defined
        
	\par Refinement of:
	- std::Callable2 <I, Vector, Vector>

	\par Associated types:
	- magnitude_type

	\par Requires:
	- std::Convertible<std::Callable2 <I, Vector, Vector>::result_type, Scalar> ;
	  result of inner product convertible to scalar to be used in expressions
	- HasConjugate < Scalar >
	- RealMagnitude < Scalar > ; the scalar value needs a real magnitude type
    */
template <typename I, typename Vector, 
          typename Scalar = typename Vector::value_type>
struct InnerProduct
  : std::Callable2<I, Vector, Vector>
{
    /// Associated type: the  real magnitude type of the scalar
    /** By default RealMagnitude<Scalar>::type */
    typename associated_type magnitude_type;
    // requires FullLessThanComparable<magnitude_type>;

    /// The arguments can be changed and the result is then the complex conjugate
    axiom ConjugateSymmetry(I inner, Vector v, Vector w)
    {
	/// inner(v, w) == conj(inner(w, v));
    }

    /// The inner product is linear in the second argument and conjugate linear in the first one
    /** The equalities are partly redundant with ConjugateSymmetry */
    axiom SequiLinearity(I inner, Scalar a, Scalar b, Vector u, Vector v, Vector w)
    {
	/// inner(v, b * w) == b * inner(v, w);

	/// inner(u, v + w) == inner(u, v) + inner(u, w);

	/// inner(a * v, w) == conj(a) * inner(v, w);

	/// inner(u + v, w) == inner(u, w) + inner(v, w);
    }

    /// The inner product of a vector with itself is not negative
    /** inner(v, v) == conj(inner(v, v)) implies inner(v, v) is representable as real */
    axiom NonNegativity(I inner, Vector v, MagnitudeType<Scalar>::type magnitude)
    {
	/// magnitude_type(inner(v, v)) >= zero(magnitude);
    }

    /// Non-degeneracy not representable with axiom
    axiom NonDegeneracy(I inner, Vector v, Vector w, Scalar s)
    {
	/// \f$\langle v, w\rangle = 0 \forall w \Leftrightarrow v = \vec{0}\f$
    }
};
#endif




#ifdef __GXX_CONCEPTS_
// A dot product is only a semantically special case of an inner product
// Questionable if we want such a concept
concept DotProduct<typename I, typename Vector, 
		   typename Scalar = typename Vector::value_type>
  : InnerProduct<I, Vector, Scalar>
{};
#else
    //! Concept DotProduct
    /*!
        Semantic requirements of dot product. The dot product is a specific inner product.

	\param I        Norm functor
	\param Vector   The the type of a vector or a collection 
        \param Scalar   The scalar over which the vector field is defined
        
	\par Refinement of:
	- InnerProduct <I, Vector, Scalar>
    */
template <typename I, typename Vector, 
	  typename Scalar = typename Vector::value_type>
struct DotProduct
  : InnerProduct<I, Vector, Scalar>
{};
#endif




// Norm induced by inner product
// Might be moved to another place later
// Definition as class and function
// Conversion from scalar to magnitude_type is covered by norm concept
template <typename I, typename Vector,
	  typename Scalar = typename Vector::value_type>
  _GLIBCXX_WHERE(InnerProduct<I, Vector, Scalar> 
		 && RealMagnitude<Scalar>)
struct induced_norm_t
{
    // Return type evtl. with macro to use concept definition
    typename magnitude_type_trait<Scalar>::type
    operator() (const I& inner, const Vector& v)
    {
	// Check whether inner product is positive real
	// assert(Scalar(abs(inner(v, v))) == inner(v, v));
	
	// Similar check while accepting small imaginary values
	// assert( (abs(inner(v, v)) - inner(v, v)) / abs(inner(v, v)) < 1e-6; )
	
	// Could also be defined with abs but that might introduce extra ops
	// typedef RealMagnitude<Scalar>::type magnitude_type;

	typedef typename magnitude_type_trait<Scalar>::type magnitude_type;
	return sqrt(static_cast<magnitude_type> (inner(v, v)));
    }
};


#if 0
template <typename I, typename Vector,
	  typename Scalar = typename Vector::value_type>
  LA_WHERE( InnerProduct<I, Vector, Scalar> 
	    && RealMagnitude<Scalar> )
magnitude_type_trait<Scalar>::type
induced_norm(const I& inner, const Vector& v)
{
    return induced_norm_t<I, Vector, Scalar>() (inner, v);
}
#endif

#ifdef __GXX_CONCEPTS__


concept HilbertSpace<typename I, typename Vector,
		     typename Scalar = typename Vector::value_type, 
		     typename N = induced_norm_t<I, Vector, Scalar> >
  : InnerProduct<I, Vector, Scalar>,
    BanachSpace<N, Vector, Scalar>
{
    axiom Consistency(Vector v)
    {
	math::induced_norm_t<I, Vector, Scalar>()(v) == N()(v);                    
    }   
};
#else
    //! Concept HilbertSpace
    /*!
        A Hilbert space is a vector space with an inner product that induces a norm

	\param I        Inner product functor
	\param Vector   The the type of a vector or a collection 
        \param Scalar   The scalar over which the vector field is defined
	\param N        Norm functor
        
	\par Refinement of:
	- InnerProduct <I, Vector, Scalar>
	- BanachSpace <N, Vector, Scalar>

	\note
	- The (expressible) requirements of Banach Space are already given in InnerProduct
	  (besides consistency of the functors).
	- A difference is that InnerProduct is not a refinement of Vectorspace
    */
template <typename I, typename Vector,
          typename Scalar = typename Vector::value_type, 
	  typename N = induced_norm_t<I, Vector, Scalar> >
struct HilbertSpace
  : InnerProduct<I, Vector, Scalar>,
    BanachSpace<N, Vector, Scalar>
{
    /// Consistency between norm and induced norm
    axiom Consistency(Vector v)
    {
	/// math::induced_norm_t<I, Vector, Scalar>()(v) == N()(v);                    
    }   
};
#endif // __GXX_CONCEPTS__

/*@}*/ // end of group Concepts

} // namespace math

#endif // LA_VECTOR_CONCEPTS_INCLUDE
