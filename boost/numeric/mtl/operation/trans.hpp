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

#ifndef MTL_TRANS_INCLUDE
#define MTL_TRANS_INCLUDE

#include <boost/mpl/if.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/algebraic_category.hpp>
#include <boost/numeric/mtl/utility/transposed_orientation.hpp>
#include <boost/numeric/mtl/utility/view_code.hpp>
#include <boost/numeric/mtl/utility/viewed_collection.hpp>
#include <boost/numeric/mtl/utility/compose_view.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/matrix/view_ref.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { 

namespace mat {

    namespace detail {

	// General case is not defined
	template <typename Value, typename AlgebraicCategory, unsigned IsConst>
	struct trans {};

	template <typename Matrix, unsigned IsConst>
	struct trans<Matrix, tag::matrix, IsConst>
	{
	    static const unsigned code= (mtl::traits::view_code<Matrix>::value | IsConst) ^ 4;
	    typedef typename mtl::traits::compose_view<code, typename mtl::traits::viewed_collection<Matrix>::type>::type result_type;
	
	    typedef typename boost::mpl::if_c<(IsConst == 1), const Matrix&, Matrix&>::type ref_type;

	    static inline result_type apply(ref_type matrix)
	    {
		return result_type(view_ref(matrix));
	    }
	};

    } // namespace detail


    template <typename Value>
    typename detail::trans<Value, typename mtl::traits::algebraic_category<Value>::type, 1>::result_type 
    inline trans(const Value& v)
    {
	vampir_trace<3041> tracer;
	return detail::trans<Value, typename mtl::traits::algebraic_category<Value>::type, 1>::apply(v);
    }

    template <typename Value>
    typename detail::trans<Value, typename mtl::traits::algebraic_category<Value>::type, 0>::result_type 
    inline trans(Value& v)
    {
	vampir_trace<3042> tracer;
	return detail::trans<Value, typename mtl::traits::algebraic_category<Value>::type, 0>::apply(v);
    }

} // namespace mtl::matrix


namespace vec {

    template <typename Vector>
    struct transposed_vector {};

    template <typename Parameters>
    struct transposed_parameters
    {
	typedef typename mtl::traits::transposed_orientation<typename Parameters::orientation>::type orientation; // switch
	typedef parameters<orientation, typename Parameters::dimension, false, typename Parameters::size_type>           type;        // not on stack!!!
    };

    template <typename Value, typename Parameters>
    struct transposed_vector<dense_vector<Value, Parameters> >
    {
	typedef dense_vector<Value, typename transposed_parameters<Parameters>::type>           type;
    };

    template <typename Value, typename Parameters>
    struct transposed_vector<strided_vector_ref<Value, Parameters> >
    {
	typedef strided_vector_ref<Value, typename transposed_parameters<Parameters>::type>     type;
    };

///Returns tranposed view of %vector v
    template <typename Vector>
    typename transposed_vector<Vector>::type const
    inline trans(const Vector& v)
    {
    vampir_trace<2037> tracer;
	typedef typename transposed_vector<Vector>::type  type;
	return type(size(v), &const_cast<Vector&>(v)[0]);
    }

    template <typename Vector>
    typename transposed_vector<Vector>::type
    inline trans(Vector& v)
    {
    vampir_trace<2038> tracer;
	typedef typename transposed_vector<Vector>::type  type;
	return type(size(v), &v[0]);
    }
}

} // mtl


#endif // MTL_TRANS_INCLUDE
