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

#ifndef MTL_CATEGORY_INCLUDE
#define MTL_CATEGORY_INCLUDE

#include <vector>

#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>

// Not elegant but necessary to treat ITL types right
#include <boost/numeric/itl/itl_fwd.hpp>


namespace mtl { namespace traits {

/// Helper for \ref category to categorize by means of \ref root
template <typename Collection> struct category_aux 
{
    typedef tag::unknown type;
};
    

/// Meta-function for categorizing MTL and external types
/** Has to be specialized for each %matrix, %vector, ...
    Extensively used for dispatching 
    @ingroup Tags
*/
template <typename Collection> struct category 
  : category_aux<Collection>
{};

// template <typename Collection> struct category 
// {
//     typedef tag::unknown type;
// };


// Const types have the same category as their non-const counterpart
template <typename T>
struct category<const T>
{
    typedef typename category<T>::type type;
};

template <typename Value, typename Parameters>
struct category<mtl::mat::dense2D<Value, Parameters> > 
{
    typedef tag::dense2D type;
};


template <typename Functor>
struct category<mtl::mat::implicit_dense<Functor> >
{
    typedef tag::implicit_dense type;
};

template <typename Value>
struct category<mtl::mat::ones_matrix<Value> >
  : public category<mtl::mat::implicit_dense<mtl::mat::ones_functor<Value> > > 
{};

template <typename Value>
struct category<mtl::mat::hilbert_matrix<Value> >
  : public category<mtl::mat::implicit_dense<mtl::mat::hilbert_functor<Value> > > 
{};

template <typename Vector1, typename Vector2>
struct category<mtl::mat::outer_product_matrix<Vector1, Vector2> >
  : public category<mtl::mat::implicit_dense<mtl::mat::outer_product_functor<Vector1, Vector2> > > 
{};


template <typename Elt, std::size_t BitMask, typename Parameters>
struct category<mtl::mat::morton_dense<Elt, BitMask, Parameters> >
{
    typedef mtl::tag::morton_dense type;
};

template <typename Elt, typename Parameters>
struct category<mtl::mat::compressed2D<Elt, Parameters> > 
{
    typedef tag::compressed2D type;
};

// should have the same tags as compressed2D
template <typename Elt, typename Parameters>
struct category<mtl::mat::coordinate2D<Elt, Parameters> > 
{
    typedef tag::compressed2D type;
};

template <typename Elt, typename Parameters>
struct category<mtl::mat::sparse_banded<Elt, Parameters> > 
{
    typedef tag::sparse_banded_matrix type;
};

template <typename Vector>
struct category<mtl::mat::multi_vector<Vector> > 
{
    typedef tag::multi_vector type;
};

template <typename Vector>
struct category<mtl::mat::multi_vector_range<Vector> > 
{
    typedef tag::multi_vector type;
};

template <typename Value>
struct category<mtl::mat::element_structure<Value> >
{
    typedef tag::element_structure type;
};

template <typename Value, typename Parameters>
struct category<mtl::mat::ell_matrix<Value, Parameters> >
{
    typedef tag::ell_matrix type;
};

template <typename T, typename Parameters>
struct category< mtl::vec::dense_vector<T, Parameters> > 
{
    typedef typename boost::mpl::if_<
	boost::is_same<typename Parameters::orientation, row_major>
      , tag::dense_row_vector 
      , tag::dense_col_vector 
    >::type type;
} ;

template <typename T, typename Parameters>
struct category< mtl::vec::strided_vector_ref<T, Parameters> > 
{
    typedef typename boost::mpl::if_<
	boost::is_same<typename Parameters::orientation, row_major>
      , tag::strided_row_vector 
      , tag::strided_col_vector 
    >::type type;
} ;

template <typename T, typename Parameters>
struct category< mtl::vec::sparse_vector<T, Parameters> > 
{
    typedef typename boost::mpl::if_<
	boost::is_same<typename Parameters::orientation, row_major>
      , tag::sparse_row_vector 
      , tag::sparse_col_vector 
    >::type type;
} ;


template <class E1, class E2, class SFunctor>
struct category< mtl::vec::vec_vec_pmop_expr<E1,E2, SFunctor> >
{
    typedef category<E1> type;
};

template <typename Functor, typename Vector> 
struct category_aux< mtl::vec::map_view<Functor, Vector> >
  : public category<Vector>
{};


// To handle std::vector in algorithms
template <typename T>
struct category< std::vector<T> >
{
    typedef tag::std_vector type;
};

namespace detail {
   
    template <typename Cat>  struct view_category       { typedef Cat                     type; };

    template <> struct view_category<tag::dense2D>      { typedef tag::dense2D_view       type; };
    template <> struct view_category<tag::morton_dense> { typedef tag::morton_view        type; };
    template <> struct view_category<tag::compressed2D> { typedef tag::compressed2D_view  type; };

    template <typename Matrix>
    struct simple_matrix_view_category
      : view_category<typename category<Matrix>::type>
    {};

} // detail


template <typename Functor, typename Matrix> 
struct category<mtl::mat::map_view<Functor, Matrix> >
  : public detail::simple_matrix_view_category<Matrix>
{};

template <typename Scaling, typename Matrix>
struct category< mtl::mat::scaled_view<Scaling, Matrix> >
    : public category< mat::map_view<tfunctor::scale<Scaling, typename Matrix::value_type>, 
					    Matrix> >
{};

// added by Hui Li
template <typename Matrix, typename RScaling>
struct category< mtl::mat::rscaled_view<Matrix,RScaling> >
    : public category< mat::map_view<tfunctor::rscale<typename Matrix::value_type,RScaling>, 
					Matrix> >
{};

// added by Hui Li
template <typename Matrix, typename Divisor>
struct category< mtl::mat::divide_by_view<Matrix,Divisor> >
    : public category< mat::map_view<tfunctor::divide_by<typename Matrix::value_type,Divisor>, 
					Matrix> >
{};

template <typename Matrix>
struct category< mtl::mat::exp_view<Matrix> >
  : category< mat::map_view<mtl::sfunctor::exp<typename Matrix::value_type>, Matrix> >
{};

// add dense to view_category
template <typename Matrix>
struct category<mtl::mat::indirect<Matrix> >
{
    typedef mtl::tag::join<typename detail::simple_matrix_view_category<Matrix>::type, mtl::tag::dense> type;
};


template <typename Matrix> 
struct category<mtl::mat::transposed_view<Matrix> >
  : public category<Matrix>
{};

// Specialize on transposed multi-vectors
template <typename Vector>
struct category< mtl::mat::transposed_view< mtl::mat::multi_vector<Vector> > >
{
    typedef tag::transposed_multi_vector type;
};

template <typename Matrix>
struct category< mat::conj_view<Matrix> >
    : public category< mat::map_view<sfunctor::conj<typename Matrix::value_type>, Matrix> >
{};

template <typename Matrix>
struct category< mat::hermitian_view<Matrix> >
  : public category< mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, 
					   mtl::mat::transposed_view<Matrix> > >
{};

// Specialize on Hermiatians of multi-vectors
template <typename Vector>
struct category< mat::hermitian_view<mtl::mat::multi_vector<Vector> > >
{
    typedef tag::hermitian_multi_vector type;
};


template <typename Matrix>
struct category< mtl::mat::banded_view<Matrix> >
  : public detail::simple_matrix_view_category<Matrix>
{};

template <typename Matrix, typename VectorIn>
struct category< mtl::vec::mat_cvec_multiplier<Matrix, VectorIn> >
{
    typedef tag::mat_cvec_multiplier type;
};

template <typename PC, typename Vector, bool Adjoint>
struct category<itl::pc::solver<PC, Vector, Adjoint> >
{
    typedef tag::unevaluated type; // might be a problem with nested matrices and vectors
};

}} // namespace mtl::traits 

#endif // MTL_CATEGORY_INCLUDE
