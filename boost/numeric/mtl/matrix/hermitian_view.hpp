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

#ifndef MTL_HERMITIAN_VIEW_INCLUDE
#define MTL_HERMITIAN_VIEW_INCLUDE

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/operation/conj.hpp>
#include <boost/numeric/mtl/operation/matrix_bracket.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>


namespace mtl { namespace mat {

template <class Matrix> 
struct hermitian_view 
  : private transposed_view<Matrix>,
    public map_view<mtl::sfunctor::conj<typename Matrix::value_type>, 
		    transposed_view<Matrix> >
{
    typedef transposed_view<Matrix>                                trans_base;
    typedef mtl::sfunctor::conj<typename Matrix::value_type>       functor_type;
    typedef map_view<functor_type, transposed_view<Matrix> >       base;
    typedef hermitian_view                                         self;
    typedef const Matrix&                                          const_ref_type;
    typedef typename Collection<Matrix>::size_type                 size_type;
    typedef typename Collection<Matrix>::value_type                value_type;

    typedef typename OrientedCollection<trans_base>::orientation   orientation; // Should not be needed because defined in Collection (bug in g++???)

    hermitian_view(const Matrix& matrix) 
      : trans_base(const_cast<Matrix&>(matrix)), 
	base(functor_type(), static_cast<trans_base&>(*this)) 
    {}
    
#if 0
    hermitian_view(boost::shared_ptr<Matrix> p)
	: trans_base(p), base(functor_type(), static_cast<trans_base&>(*this))
    {}
#endif	

    typename base::value_type operator()(size_type r, size_type c) const { return base::operator()(r, c); }

    operations::bracket_proxy<self, const self&, value_type> 
    operator[] (size_type r) const
    {
	return operations::bracket_proxy<self, const self&, value_type>(*this, r);
    }

    friend size_type inline num_rows(const self& A) { return num_rows((const base&)(A)); }
    friend size_type inline num_cols(const self& A) { return num_cols((const base&)(A)); }
 
    const_ref_type const_ref() const 
    { 
	// make two statements because nvcc cannot handle ref.ref
	const transposed_view<Matrix>& r1= base::ref;
	return r1.ref; 
    }

    size_type nnz() const { return base::nnz(); }

    friend inline std::ostream& operator<<(std::ostream& os, const self& A) { return os << (const base&)(A); }
};

// If not defined ambigous between map_view and transposed_view
template <class Matrix> 
inline std::size_t size(const hermitian_view<Matrix>& A)
{  
    return num_rows(A) * num_rows(A); 
}

// TBD submatrix of Hermitian (not trivial)


}} // namespace mtl::matrix



// Traits for Hermitian views
namespace mtl { namespace traits {

template <typename Matrix>
struct row< mtl::mat::hermitian_view<Matrix> >
  : public row< mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, 
				      mtl::mat::transposed_view<Matrix> > >
{};

template <typename Matrix>
struct col< mtl::mat::hermitian_view<Matrix> >
    : public col< mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, 
				   mtl::mat::transposed_view<Matrix> > >
{};

template <typename Matrix>
struct const_value< mtl::mat::hermitian_view<Matrix> >
    : public const_value< mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, 
					   mtl::mat::transposed_view<Matrix> > >
{};

template <typename Tag, typename Matrix>
struct range_generator< Tag, mtl::mat::hermitian_view<Matrix> >
    : public range_generator< Tag, mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, 
						    mtl::mat::transposed_view<Matrix> > >
{};

template <typename Matrix>
struct range_generator< tag::major, mtl::mat::hermitian_view<Matrix> >
    : public range_generator< tag::major, mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, 
							   mtl::mat::transposed_view<Matrix> > >
{};


}} // mtl::traits

#endif // MTL_HERMITIAN_VIEW_INCLUDE
