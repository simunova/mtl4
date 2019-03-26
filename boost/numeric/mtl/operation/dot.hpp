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

#ifndef MTL_DOT_INCLUDE
#define MTL_DOT_INCLUDE


#include <boost/numeric/mtl/concept/std_concept.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/conj.hpp>
#include <boost/numeric/meta_math/loop1.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/mtl/utility/omp_size_type.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>

namespace mtl { 

    namespace vec {

	namespace detail {
	    
	    // Result type of dot product
	    template <typename Vector1, typename Vector2>
	    struct dot_result
	    {
		typedef typename Multiplicable<typename Collection<Vector1>::value_type,
					       typename Collection<Vector2>::value_type>::result_type type;
	    };

	    // Whether or not conjugating first argument
	    struct without_conj
	    {
		template <typename Value>
		Value operator()(const Value& v) { return v; }
	    };

	    struct with_conj
	    {
		template <typename Value>
		typename mtl::sfunctor::conj<Value>::result_type
		operator() (const Value& v)
		{
		    using mtl::conj;
		    return conj(v);
		}
	    };
	}

	namespace sfunctor {
	    
	    template <unsigned long Index0, unsigned long Max0>
	    struct dot_aux
		: public meta_math::loop1<Index0, Max0>
	    {
		typedef meta_math::loop1<Index0, Max0>                                    base;
		typedef dot_aux<base::next_index0, Max0>                                  next_t;
		
		template <typename Value, typename Vector1, typename Vector2, typename Size, typename ConjOpt>
		static inline void 
		apply(Value& tmp00, Value& tmp01, Value& tmp02, Value& tmp03, Value& tmp04, 
		      Value& tmp05, Value& tmp06, Value& tmp07, 
		      const Vector1& v1, const Vector2& v2, Size i, ConjOpt conj_opt)
		{
		    // vampir_trace<9901> tracer;
		    tmp00+= conj_opt(v1[ i + base::index0 ]) * v2[ i + base::index0 ];
		    next_t::apply(tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp00,
				  v1, v2, i, conj_opt);
		}
	    };


	    template <unsigned long Max0>
	    struct dot_aux<Max0, Max0>
	    {
		typedef meta_math::loop1<Max0, Max0>                                      base;
		
		template <typename Value, typename Vector1, typename Vector2, typename Size, typename ConjOpt>
		static inline void 
		apply(Value& tmp00, Value&, Value&, Value&, Value&, Value&, Value&, Value&, 
		      const Vector1& v1, const Vector2& v2, Size i, ConjOpt conj_opt)
		{
		    tmp00+= conj_opt(v1[ i + base::index0 ]) * v2[ i + base::index0 ];
		}
	    };


	    template <unsigned long Unroll>
	    struct dot
	    {
		template <typename Vector1, typename Vector2, typename ConjOpt>
		typename detail::dot_result<Vector1, Vector2>::type
		static inline apply(const Vector1& v1, const Vector2& v2, ConjOpt conj_opt)
		{
		    MTL_STATIC_ASSERT((Unroll >= 1), "Unroll size must be at least 1.");
		    // MTL_STATIC_ASSERT((Unroll <= 8), "Maximal unrolling is 8."); // Might be relaxed in future versions

		    vampir_trace<2003> tracer;
		    MTL_THROW_IF(mtl::size(v1) != mtl::size(v2), incompatible_size());
		    typedef typename detail::dot_result<Vector1, Vector2>::type  value_type;
		    		    
#                 ifdef MTL_WITH_OPENMP 
		    value_type dummy, z= math::zero(dummy), result= z;
		    typedef typename mtl::traits::omp_size_type<typename Collection<Vector1>::size_type>::type size_type;
		    size_type  i_max= mtl::size(v1), i_block= Unroll * (i_max / Unroll);


                    #pragma omp parallel
		    {

			vampir_trace<8001> tracer;
			value_type tmp00= z, tmp01= z, tmp02= z, tmp03= z, tmp04= z, tmp05= z, tmp06= z, tmp07= z;

			#pragma omp for
			for (size_type i= 0; i < i_block; i+= Unroll)
			    dot_aux<1, Unroll>::apply(tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, v1, v2, i, conj_opt);

			#pragma omp critical
			    result+= ((tmp00 + tmp01) + (tmp02 + tmp03)) + ((tmp04 + tmp05) + (tmp06 + tmp07));
		    }
		    for (size_type i= i_block; i < i_max; i++) 
			result+= conj_opt(v1[i]) * v2[i];

		    return result;
#                 else
		    typedef typename Collection<Vector1>::size_type              size_type;

		    value_type dummy, z= math::zero(dummy), tmp00= z, tmp01= z, tmp02= z, tmp03= z, tmp04= z,
			       tmp05= z, tmp06= z, tmp07= z;
		    size_type  i_max= mtl::size(v1), i_block= Unroll * (i_max / Unroll);
		    
		    for (size_type i= 0; i < i_block; i+= Unroll)
			dot_aux<1, Unroll>::apply(tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, v1, v2, i, conj_opt);
		    
		    for (size_type i= i_block; i < i_max; i++) 
			tmp00+= conj_opt(v1[i]) * v2[i];
		    return ((tmp00 + tmp01) + (tmp02 + tmp03)) + ((tmp04 + tmp05) + (tmp06 + tmp07));
#                 endif
		}


	    };
	}

	template <typename Vector1, typename Vector2, typename ConjOpt>
	typename detail::dot_result<Vector1, Vector2>::type
	inline dot_simple(const Vector1& v1, const Vector2& v2, ConjOpt conj_opt)
	{
	    vampir_trace<2040> tracer;
	    typedef typename Collection<Vector1>::size_type              size_type;
	    typedef typename detail::dot_result<Vector1, Vector2>::type  value_type;

	    value_type dummy, s= math::zero(dummy);
	    for (size_type i= 0, i_max= mtl::size(v1); i < i_max; ++i)
		s+= conj_opt(v1[i]) * v2[i];
	    return s;
	}

	template <unsigned long Unroll, typename Vector1, typename Vector2, typename ConjOpt>
	struct dot_class
	{
	    typedef typename detail::dot_result<Vector1, Vector2>::type result_type;
	    dot_class(const Vector1& v1, const Vector2& v2) : v1(v1), v2(v2) {}

	    operator result_type() const { return sfunctor::dot<Unroll>::apply(v1, v2, ConjOpt()); }
	    
	    const Vector1& v1;
	    const Vector2& v2;
	};

	template <typename Vector1, typename Vector2, typename ConjOpt>
	struct dot_class<1, Vector1, Vector2, ConjOpt>
	{
	    typedef typename detail::dot_result<Vector1, Vector2>::type result_type;
	    dot_class(const Vector1& v1, const Vector2& v2) : v1(v1), v2(v2) {}

	    operator result_type() const { return dot_simple(v1, v2, ConjOpt()); }
	    
	    const Vector1& v1;
	    const Vector2& v2;
	};
	
	/// Lazy dot product
	/** It is automatically evaluated when (implicitly) converted to result_type which doesn't work in template expressions.
	    Can be used for source-to-source transformations. **/
	template <typename Vector1, typename Vector2>
	dot_class<4, Vector1, Vector2, detail::with_conj>
	inline lazy_dot(const Vector1& v1, const Vector2& v2)
	{
	    return dot_class<4, Vector1, Vector2, detail::with_conj>(v1, v2);
	}

	template <unsigned long Unroll, typename Vector1, typename Vector2>
	dot_class<Unroll, Vector1, Vector2, detail::with_conj>
	inline lazy_dot(const Vector1& v1, const Vector2& v2)
	{
	    return dot_class<Unroll, Vector1, Vector2, detail::with_conj>(v1, v2);
	}

	template <typename Vector1, typename Vector2>
	dot_class<4, Vector1, Vector2, detail::without_conj>
	inline lazy_dot_real(const Vector1& v1, const Vector2& v2)
	{
	    return dot_class<4, Vector1, Vector2, detail::without_conj>(v1, v2);
	}

	template <unsigned long Unroll, typename Vector1, typename Vector2>
	dot_class<Unroll, Vector1, Vector2, detail::without_conj>
	inline lazy_dot_real(const Vector1& v1, const Vector2& v2)
	{
	    return dot_class<Unroll, Vector1, Vector2, detail::without_conj>(v1, v2);
	}

	/// Dot product defined as hermitian(v) * w
	/** Unrolled four times by default **/
	template <typename Vector1, typename Vector2>
	typename detail::dot_result<Vector1, Vector2>::type
	inline dot(const Vector1& v1, const Vector2& v2)
	{
	    // return dot_simple(v1, v2, detail::with_conj());
	    return sfunctor::dot<4>::apply(v1, v2, detail::with_conj());
	}

	/// Dot product with user-specified unrolling defined as hermitian(v) * w
	template <unsigned long Unroll, typename Vector1, typename Vector2>
	typename detail::dot_result<Vector1, Vector2>::type
	inline dot(const Vector1& v1, const Vector2& v2)
	{
	    return sfunctor::dot<Unroll>::apply(v1, v2, detail::with_conj());
	}
	/// Dot product without conjugate defined as trans(v) * w
	/** Unrolled four times by default **/
	template <typename Vector1, typename Vector2>
	typename detail::dot_result<Vector1, Vector2>::type
	inline dot_real(const Vector1& v1, const Vector2& v2)
	{
	    return sfunctor::dot<4>::apply(v1, v2, detail::without_conj());
	}

	/// Dot product without conjugate with user-specified unrolling defined as trans(v) * w
	template <unsigned long Unroll, typename Vector1, typename Vector2>
	typename detail::dot_result<Vector1, Vector2>::type
	inline dot_real(const Vector1& v1, const Vector2& v2)
	{
	    return sfunctor::dot<Unroll>::apply(v1, v2, detail::without_conj());
	}


    } // namespace vector
    
    using vec::dot;
    using vec::dot_real;
    using vec::lazy_dot;
    using vec::lazy_dot_real;

} // namespace mtl

#endif // MTL_DOT_INCLUDE
