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

#ifndef MTL_MATRIX_IMPLICIT_DENSE_INCLUDE
#define MTL_MATRIX_IMPLICIT_DENSE_INCLUDE

#include <vector>
#include <boost/numeric/linear_algebra/inverse.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/matrix/crtp_base_matrix.hpp>
#include <boost/numeric/mtl/matrix/mat_expr.hpp>
#include <boost/numeric/mtl/concept/std_concept.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>

#ifdef MTL_HAS_MPI
#   include <boost/numeric/mtl/matrix/distributed.hpp>
#   include <boost/numeric/mtl/matrix/inserter.hpp>
#   include <boost/mpi/collectives.hpp>
#endif

namespace mtl { namespace mat {

template <typename Functor>
class implicit_dense
  : public const_crtp_base_matrix< implicit_dense<Functor>,
				   typename Functor::result_type, typename Functor::size_type >,
    public mat_expr< implicit_dense<Functor> >
{
    typedef implicit_dense                             self;
  public:
    typedef mtl::tag::row_major                        orientation; // for completeness
    typedef typename Functor::result_type              value_type;
    typedef typename Functor::result_type              const_reference;
    typedef typename Functor::size_type                size_type;

    // Should not be needed, to be removed after transposed_view is cleaned up
    typedef index::c_index                             index_type;
    typedef mtl::traits::detail::matrix_element_key<self> key_type;
    typedef mtl::non_fixed::dimensions                 dim_type;

    explicit implicit_dense (const Functor& functor) : my_functor(functor) {}

    value_type operator() (size_type r, size_type c) const { return my_functor(r, c); }

    size_type nnz() const { return dim1() * dim2(); }
    size_type dim1() const { return num_rows(my_functor); }
    size_type dim2() const { return num_cols(my_functor); }

    friend size_type inline num_rows(const self& A) { return num_rows(A.my_functor); }
    friend size_type inline num_cols(const self& A) { return num_cols(A.my_functor); }

    Functor& functor() { return my_functor; }
    Functor const& functor() const { return my_functor; }

  private:
    Functor           my_functor;    
};

template <typename Functor>
typename Functor::size_type inline size(const implicit_dense<Functor>& A) 
{ 
    return num_rows(A) * num_cols(A);
}

// ==========
// Sub matrix
// ==========

// To do later

}} // namespace mtl::matrix


namespace mtl { namespace traits {

    // ================
    // Range generators
    // For cursors
    // ================

    template <typename Functor>
    struct range_generator<glas::tag::row, mtl::mat::implicit_dense<Functor> >
      : detail::all_rows_range_generator<mtl::mat::implicit_dense<Functor>, complexity_classes::linear>
    {};

    template <typename Functor>
    struct range_generator<glas::tag::major, mtl::mat::implicit_dense<Functor> >
      : range_generator<glas::tag::row, mtl::mat::implicit_dense<Functor> >
    {};

    template <typename Functor>
    struct range_generator<glas::tag::nz,
			   detail::sub_matrix_cursor<mtl::mat::implicit_dense<Functor>, glas::tag::row, 2> > 
      : detail::all_cols_in_row_range_generator<detail::sub_matrix_cursor<mtl::mat::implicit_dense<Functor>, glas::tag::row, 2> >
    {};

    template <typename Functor>
    struct range_generator<glas::tag::col, mtl::mat::implicit_dense<Functor> >
      : detail::all_cols_range_generator<mtl::mat::implicit_dense<Functor>, complexity_classes::linear>
    {};

    template <typename Functor>
    struct range_generator<glas::tag::nz,
			   detail::sub_matrix_cursor<mtl::mat::implicit_dense<Functor>, glas::tag::col, 2> > 
      : detail::all_rows_in_col_range_generator<detail::sub_matrix_cursor<mtl::mat::implicit_dense<Functor>, glas::tag::col, 2> >
    {};


    template <typename Tag, typename Value>
    struct range_generator<Tag, mtl::mat::ones_matrix<Value> >
      : public range_generator<Tag, mtl::mat::implicit_dense<mtl::mat::ones_functor<Value> > > 
    {};

    template <typename Value>
    struct range_generator<glas::tag::major, mtl::mat::ones_matrix<Value> >
      : public range_generator<glas::tag::major, mtl::mat::implicit_dense<mtl::mat::ones_functor<Value> > > 
    {};

    template <typename Tag, typename Value>
    struct range_generator<Tag, mtl::mat::hilbert_matrix<Value> >
      : public range_generator<Tag, mtl::mat::implicit_dense<mtl::mat::hilbert_functor<Value> > > 
    {};

    template <typename Value>
    struct range_generator<glas::tag::major, mtl::mat::hilbert_matrix<Value> >
      : public range_generator<glas::tag::major, mtl::mat::implicit_dense<mtl::mat::hilbert_functor<Value> > > 
    {};

    template <typename Tag, typename Vector1, typename Vector2>
    struct range_generator<Tag, mtl::mat::outer_product_matrix<Vector1, Vector2> >
      : public range_generator<Tag, mtl::mat::implicit_dense<mtl::mat::outer_product_functor<Vector1, Vector2> > > 
    {};

    template <typename Vector1, typename Vector2>
    struct range_generator<glas::tag::major, mtl::mat::outer_product_matrix<Vector1, Vector2> >
      : public range_generator<glas::tag::major, mtl::mat::implicit_dense<mtl::mat::outer_product_functor<Vector1, Vector2> > > 
    {};

}} // mtl::traits

// =============
// Some functors
// =============


namespace mtl { namespace mat {

template <typename Value= int>
class ones_functor
{
    typedef ones_functor    self;
  public:
    typedef std::size_t    size_type;
    typedef Value          result_type;

    ones_functor(size_type nr, size_type nc) : nr(nr), nc(nc) {}

    friend size_type inline num_rows(const self& A) { return A.nr; }
    friend size_type inline num_cols(const self& A) { return A.nc; }

    result_type operator()(size_type, size_type) const { return Value(1); }

  private:
    size_type nr, nc;
};

    
template <typename Value= double>
class hilbert_functor
{
    typedef hilbert_functor    self;
  public:
    typedef std::size_t    size_type;
    typedef Value          result_type;

    hilbert_functor(size_type nr, size_type nc) : nr(nr), nc(nc) {}

    friend size_type inline num_rows(const self& A) { return A.nr; }
    friend size_type inline num_cols(const self& A) { return A.nc; }

    result_type operator()(size_type r, size_type c) const 
    { 
	using math::reciprocal;  
	return reciprocal(Value(r + c + 1)); 
    }
  private:
    size_type nr, nc;
};

    
template <typename Vector1, typename Vector2>
class outer_product_functor
{
    typedef outer_product_functor    self;
  public:
    typedef std::size_t              size_type;
    typedef typename Multiplicable<typename Collection<Vector1>::value_type,
				   typename Collection<Vector2>::value_type>::result_type   result_type;

    outer_product_functor(const Vector1& v1, const Vector2& v2) : my_v1(v1), my_v2(v2) {}
    outer_product_functor(size_type r, size_type c) : my_v1(r), my_v2(c) {}

    friend size_type inline num_rows(const self& A) { return mtl::size(A.my_v1); }
    friend size_type inline num_cols(const self& A) { return mtl::size(A.my_v2); }

    result_type operator()(size_type r, size_type c) const { return my_v1[r] * my_v2[c]; }

    Vector1&       v1()       { return my_v1; }
    Vector1 const& v1() const { return my_v1; }
    Vector2&       v2()       { return my_v2; }
    Vector2 const& v2() const { return my_v2; }

  private:
    Vector1  my_v1; // keeps copy
    Vector2  my_v2;
};

// ======================
// Some implicit matrices
// ======================


template <typename Value= int>
class ones_matrix
  : public implicit_dense<ones_functor<Value> >
{
  public:
    typedef ones_functor<Value>                functor_type;
    typedef typename functor_type::size_type  size_type;
    typedef implicit_dense<functor_type>      base;

    ones_matrix(size_type r, size_type c) : base(functor_type(r, c)) {}
};


template <typename Value= double>
class hilbert_matrix
  : public implicit_dense<hilbert_functor<Value> >
{
  public:
    typedef hilbert_functor<Value>            functor_type;
    typedef typename functor_type::size_type  size_type;
    typedef implicit_dense<functor_type>      base;

    hilbert_matrix(size_type r, size_type c) : base(functor_type(r, c)) {}
};


template <typename Vector1, typename Vector2>
class outer_product_matrix
  : public implicit_dense<outer_product_functor<Vector1, Vector2> >
{
    typedef outer_product_matrix self;
  public:
    typedef outer_product_functor<Vector1, Vector2>  functor_type;
    typedef typename functor_type::size_type  size_type;
    typedef implicit_dense<functor_type>      base;

    outer_product_matrix(const Vector1& v1, const Vector2& v2) : base(functor_type(v1, v2)) {}
    outer_product_matrix(size_type r, size_type c) : base(functor_type(r, c)) {}

    Vector1&       v1()       { return this->functor().v1(); }
    Vector1 const& v1() const { return this->functor().v1(); }
    Vector2&       v2()       { return this->functor().v2(); }
    Vector2 const& v2() const { return this->functor().v2(); }
};

}} // namespace mtl::matrix

#endif // MTL_MATRIX_IMPLICIT_DENSE_INCLUDE
