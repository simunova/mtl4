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

#ifndef MTL_TRAITS_EVAL_INCLUDE
#define MTL_TRAITS_EVAL_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/matrix/crtp_base_matrix.hpp>
#include <boost/numeric/mtl/matrix/mat_expr.hpp>

namespace mtl { namespace traits {


template <typename T>
struct eval {};

#if 0 // To be done later
template <typename Value, typename Parameter>
struct eval< mtl::dense_vector<Value, Parameter> >
{};


template <typename Value1, typename Vector>
struct eval< mtl::scaled_view<Value1, Vector> > 
{};

template <typename Value1, typename Vector>
struct eval< mtl::rscaled_view<Value1, Vector> > 
{};
#endif


    namespace impl {
	
	template<typename T>
	struct eval_self_ref 
	{
	    typedef const T&     const_reference;

	    explicit eval_self_ref(const T& ref) : ref(ref) {}

	    const_reference value() const { return ref; }

	    const T& ref;
	};
    }


template <typename Value, typename Parameter>
struct eval< mtl::mat::dense2D<Value, Parameter> >
    : public impl::eval_self_ref< mtl::mat::dense2D<Value, Parameter> >
{
    eval(const mtl::mat::dense2D<Value, Parameter>& ref)
	: impl::eval_self_ref< mtl::mat::dense2D<Value, Parameter> >(ref)
    {}
};

template <typename Value, std::size_t Mask, typename Parameter>
struct eval< mtl::mat::morton_dense<Value, Mask, Parameter> >
    : public impl::eval_self_ref< mtl::mat::morton_dense<Value, Mask, Parameter> >
{
    eval(const mtl::mat::morton_dense<Value, Mask, Parameter>& ref)
	: impl::eval_self_ref< mtl::mat::morton_dense<Value, Mask, Parameter> >(ref)
    {}
};

template <typename Value, typename Parameter>
struct eval< mtl::mat::compressed2D<Value, Parameter> >
    : public impl::eval_self_ref< mtl::mat::compressed2D<Value, Parameter> >
{
    eval(const mtl::mat::compressed2D<Value, Parameter>& ref)
	: impl::eval_self_ref< mtl::mat::compressed2D<Value, Parameter> >(ref)
    {}
};




#if 0 // only dummy
template <typename E1, typename E2>
struct eval< mtl::mat::mat_mat_asgn_expr<E1, E2> > 
{};
#endif


template <typename E1, typename E2>
struct eval< mtl::mat::mat_mat_plus_expr<E1, E2> > 
    : public impl::eval_self_ref< mtl::mat::mat_mat_plus_expr<E1, E2> > 
{
    eval(const mtl::mat::mat_mat_plus_expr<E1, E2>& ref)
	: impl::eval_self_ref< mtl::mat::mat_mat_plus_expr<E1, E2> >(ref)
    {}
};

template <typename E1, typename E2>
struct eval< mtl::mat::mat_mat_minus_expr<E1, E2> > 
    : public impl::eval_self_ref< mtl::mat::mat_mat_minus_expr<E1, E2> > 
{
    eval(const mtl::mat::mat_mat_minus_expr<E1, E2>& ref)
	: impl::eval_self_ref< mtl::mat::mat_mat_minus_expr<E1, E2> >(ref)
    {}
};

template <typename E1, typename E2>
struct eval< mtl::mat::mat_mat_ele_times_expr<E1, E2> > 
    : public impl::eval_self_ref< mtl::mat::mat_mat_ele_times_expr<E1, E2> > 
{
    eval(const mtl::mat::mat_mat_ele_times_expr<E1, E2>& ref)
	: impl::eval_self_ref< mtl::mat::mat_mat_ele_times_expr<E1, E2> >(ref)
    {}
};


template <typename E1, typename E2>
struct eval< mtl::mat::mat_mat_times_expr<E1, E2> > 
{
    // Needs dramatic improvement!!! Only for testing!!!
    typedef mat::dense2D<double>  matrix_type;
    typedef const matrix_type&       const_reference;


    explicit eval(const mtl::mat::mat_mat_times_expr<E1, E2>& expr)
    	: prod(expr.first * expr.second)
    {}

    const_reference value() { return prod; }

private:
    matrix_type  prod;
};




template <typename Value1, typename Matrix>
struct eval< mtl::mat::scaled_view<Value1, Matrix> > 
{};

template <typename Value1, typename Matrix>
struct eval< mtl::mat::rscaled_view<Value1, Matrix> > 
{};


template <typename T>
eval<T> inline evaluate(const T& ref)
{
    return eval<T>(ref);
}


} // namespace traits

namespace mat {
    using mtl::traits::evaluate;
}

namespace vec {
    using mtl::traits::evaluate;
}


} // namespace mtl

#endif // MTL_TRAITS_EVAL_INCLUDE
