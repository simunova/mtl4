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

#ifndef MTL_COPY_INCLUDE
#define MTL_COPY_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/detail/index.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/flatcat.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/utility/updater_to_assigner.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/operation/update.hpp>
#include <boost/numeric/mtl/operation/print.hpp>
#include <boost/numeric/mtl/operation/crop.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <iostream>
#include <limits>

namespace mtl {
	
    namespace detail {

	// Set Destination matrix to zero when source is sparse 
	// (otherwise everything is overwritten anyway)
	template <typename MatrixDest>
	inline void zero_with_sparse_src(MatrixDest& dest, tag::flat<tag::sparse>)
	{
	    set_to_zero(dest);
	}
	
	template <typename MatrixDest>
	inline void zero_with_sparse_src(MatrixDest&, tag::universe) {} 

	// Adapt inserter size to operation
	template <typename Updater> struct copy_inserter_size {};
	
	// Specialization for store
	template <typename Value>
	struct copy_inserter_size< operations::update_store<Value> >
	{
	    template <typename MatrixSrc, typename MatrixDest>
	    static inline int apply(const MatrixSrc& src, const MatrixDest& dest)
	    {
		// std::cout << "nnz = " << src.nnz() << ", dim1 = " << dest.dim1() << "\n";
		return int(src.nnz() * 1.2 / dest.dim1());
	    }
	};

	struct sum_of_sizes
	{
	    template <typename MatrixSrc, typename MatrixDest>
	    static inline int apply(const MatrixSrc& src, const MatrixDest& dest)
	    {	return int((src.nnz() + dest.nnz()) * 1.2 / dest.dim1()); }
	};
	    	
	// Specialization for plus and minus
	template <typename Value> struct copy_inserter_size< operations::update_plus<Value> > : sum_of_sizes {};
	template <typename Value> struct copy_inserter_size< operations::update_minus<Value> > : sum_of_sizes {};

    } // namespace detail


    template <typename Updater, typename MatrixSrc, typename MatrixDest>
    inline void gen_matrix_copy(const MatrixSrc& src, MatrixDest& dest, bool with_reset)
    {
	vampir_trace<3002> tracer;
	MTL_THROW_IF(num_rows(src) != num_rows(dest) || num_cols(src) != num_cols(dest), incompatible_size());

	if (with_reset)
	    detail::zero_with_sparse_src(dest, traits::sparsity_flatcat<MatrixSrc>()); 
	
	typename traits::row<MatrixSrc>::type             row(src); 
	typename traits::col<MatrixSrc>::type             col(src); 
	typename traits::const_value<MatrixSrc>::type     value(src); 
	typedef typename traits::range_generator<tag::major, MatrixSrc>::type  cursor_type;

	// std::cout << "Slot size is " << detail::copy_inserter_size<Updater>::apply(src, dest) << "\n";
	mat::inserter<MatrixDest, Updater>   ins(dest, detail::copy_inserter_size<Updater>::apply(src, dest));
	for (cursor_type cursor = mtl::begin<tag::major>(src), cend = mtl::end<tag::major>(src); 
	     cursor != cend; ++cursor) {
	    // std::cout << dest << '\n';
	    
	    typedef typename traits::range_generator<tag::nz, cursor_type>::type icursor_type;
	    for (icursor_type icursor = mtl::begin<tag::nz>(cursor), icend = mtl::end<tag::nz>(cursor); 
		 icursor != icend; ++icursor) {
		//std::cout << "in " << row(*icursor) << ", " << col(*icursor) << " insert " << value(*icursor) << '\n';
		ins(row(*icursor), col(*icursor)) << value(*icursor); }
	}
    }

    // Specialization for multi_vector
    template <typename Updater, typename MatrixSrc, typename Vector>
    inline void gen_matrix_copy(const MatrixSrc& src, mtl::mat::multi_vector<Vector>& dest, bool)
    {
	MTL_THROW_IF(num_rows(src) != num_rows(dest) || num_cols(src) != num_cols(dest), incompatible_size());
	typedef typename mtl::traits::updater_to_assigner<Updater>::type Assigner;

	for (std::size_t i= 0, n= num_cols(src); i < n; ++i)
	    Assigner::first_update(dest.vector(i), src.vector(i));
    }

    namespace {
#    ifdef __clang__
#     pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#     pragma clang diagnostic ignored "-Wunused-function"
#    endif
		template <typename T>
		inline T inc_wo_over(T i) 
		{ return i == std::numeric_limits<T>::max() ? i : i+1; }

		template <typename T>
		inline T negate_wo_over(T i) 
		{ return i == std::numeric_limits<T>::min() ? std::numeric_limits<T>::max() : -i; }
    }

    template <typename Updater, typename ValueSrc, typename Para, typename ValueDest>
    typename boost::enable_if<boost::is_same<Updater, operations::update_store<ValueDest> > >::type
    inline gen_matrix_copy(const mat::banded_view<mtl::mat::compressed2D<ValueSrc, Para> >& src, mtl::mat::compressed2D<ValueDest, Para>& dest, bool)
    {
	vampir_trace<3061> tracer;
	typedef typename Para::size_type size_type;
	dest.change_dim(num_rows(src), num_cols(src)); // contains make_empty
	set_to_zero(dest);
	const mtl::mat::compressed2D<ValueSrc, Para>  &sref= src.ref;
	const std::vector<size_type>        &sstarts= sref.ref_major(), &sindices= sref.ref_minor();
	long first, last;
	if (traits::is_row_major<Para>::value) {
	    first= src.get_begin();
	    last= src.get_end();
	} else {
	    first= inc_wo_over(negate_wo_over(src.get_end()));
	    last=  inc_wo_over(negate_wo_over(src.get_begin()));
	}

	long jd= 0, j_end= long(sstarts[0]);
	for (long i= 0, i_end= long(src.dim1()), f= first, l= last; i < i_end; ++i) {
	    dest.ref_major()[i]= jd;
	    long j= j_end;
	    j_end= long(sstarts[i+1]);
	    while (j < j_end && long(sindices[j]) < f) j++;
	    while (j < j_end && long(sindices[j]) < l) jd++, j++;
	    f= inc_wo_over(f);
	    l= inc_wo_over(l);
	}
	dest.ref_major()[src.dim1()]= jd;
	dest.set_nnz(jd); // resizes indices and data

	for (long i= 0, i_end= long(src.dim1()), jd= 0, j_end= long(sstarts[0]); i < i_end; ++i) {
	    dest.ref_major()[i]= jd;
	    long j= j_end;
	    j_end= long(sstarts[i+1]);
	    while (j < j_end && long(sindices[j]) < first) j++;
	    while (j < j_end && long(sindices[j]) < last) {
		dest.ref_minor()[jd]= sindices[j];
		dest.data[jd++]= sref.data[j++];
	    }
	    first= inc_wo_over(first);
	    last= inc_wo_over(last);	    
	}
    }

	    
    /// Copy matrix \p src into matrix \p dest
    template <typename MatrixSrc, typename MatrixDest>
    inline void matrix_copy(const MatrixSrc& src, MatrixDest& dest)
    {
	gen_matrix_copy< operations::update_store<typename MatrixDest::value_type> >(src, dest, true);
    }
    

    /// Add matrix \p src to matrix \p dest in copy-like style
    template <typename MatrixSrc, typename MatrixDest>
    inline void matrix_copy_plus(const MatrixSrc& src, MatrixDest& dest)
    {
	gen_matrix_copy< operations::update_plus<typename MatrixDest::value_type> >(src, dest, false);
    }
	
    /// Subtract matrix \p src from matrix \p dest in copy-like style
    template <typename MatrixSrc, typename MatrixDest>
    inline void matrix_copy_minus(const MatrixSrc& src, MatrixDest& dest)
    {
	gen_matrix_copy< operations::update_minus<typename MatrixDest::value_type> >(src, dest, false);
    }
	
    /// Multiply matrix \p src element-wise with matrix \p dest in copy-like style
    template <typename MatrixSrc, typename MatrixDest>
    inline void matrix_copy_ele_times(const MatrixSrc& src, MatrixDest& dest)
    {
	vampir_trace<3001> tracer;
	MTL_THROW_IF(num_rows(src) != num_rows(dest) || num_cols(src) != num_cols(dest), incompatible_size());

	typename traits::row<MatrixDest>::type             row(dest); 
	typename traits::col<MatrixDest>::type             col(dest); 
	typename traits::value<MatrixDest>::type           value(dest); 
	typedef typename traits::range_generator<tag::major, MatrixDest>::type  cursor_type;
	typedef typename traits::range_generator<tag::nz, cursor_type>::type icursor_type;
	
	for (cursor_type cursor = begin<tag::major>(dest), cend = end<tag::major>(dest); cursor != cend; ++cursor)
	    for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); icursor != icend; ++icursor)
		value(*icursor, value(*icursor) * src[row(*icursor)][col(*icursor)]);
#if 0   // copy would result in a*0 = a and 0*b = b!!!!
	gen_matrix_copy< operations::update_times<typename MatrixDest::value_type> >(src, dest, false);
#endif
	crop(dest);
    }

       
    template <typename MatrixSrc, typename MatrixDest>
    inline void copy(const MatrixSrc& src, tag::flat<tag::matrix>, MatrixDest& dest, tag::flat<tag::matrix>)
	// inline void copy(const MatrixSrc& src, tag::matrix_expr, MatrixDest& dest, tag::matrix)
    {
	return matrix_copy(src, dest);
    }

  

    template <typename Updater, typename VectorSrc, typename VectorDest>
    inline void gen_vector_copy(const VectorSrc& src, VectorDest& dest, bool with_reset)
    {
	// Works only with dense vectors as dest !!!!! (source could be sparse)
	// Needs vector inserter
	vampir_trace<2001> tracer;

	MTL_STATIC_ASSERT((boost::is_same<typename ashape::ashape<VectorSrc>::type,
					  typename ashape::ashape<VectorDest>::type>::value), "Source and target must have the same algebraic shape.");

	MTL_THROW_IF(size(src) != size(dest), incompatible_size());

	if (with_reset)
	    detail::zero_with_sparse_src(dest, typename traits::category<VectorSrc>::type());
	
	typename traits::index<VectorSrc>::type           index(src); 
	typename traits::const_value<VectorSrc>::type     value(src); 

	typedef typename traits::range_generator<tag::nz, VectorSrc>::type  cursor_type;
	for (cursor_type cursor = begin<tag::nz>(src), cend = end<tag::nz>(src); 
	     cursor != cend; ++cursor)
	    Updater()(dest[index(*cursor)], value(*cursor));
    }
	    
    /// Copy vector \p src into vector \p dest
    template <typename VectorSrc, typename VectorDest>
    inline void vector_copy(const VectorSrc& src, VectorDest& dest)
    {
	gen_vector_copy< operations::update_store<typename VectorDest::value_type> >(src, dest, true);
    }
    

    /// Add vector \p src to vector \p dest in copy-like style
    template <typename VectorSrc, typename VectorDest>
    inline void vector_copy_plus(const VectorSrc& src, VectorDest& dest)
    {
	gen_vector_copy< operations::update_plus<typename VectorDest::value_type> >(src, dest, false);
    }
	
    /// Subtract vector \p src from vector \p dest in copy-like style
    template <typename VectorSrc, typename VectorDest>
    inline void vector_copy_minus(const VectorSrc& src, VectorDest& dest)
    {
	gen_vector_copy< operations::update_minus<typename VectorDest::value_type> >(src, dest, false);
    }
	

       
    template <typename VectorSrc, typename VectorDest>
    inline void copy(const VectorSrc& src, tag::flat<tag::vector>, VectorDest& dest, tag::flat<tag::vector>)	
    {
	return vector_copy(src, dest);
    }


    template <typename CollSrc, typename CollDest>
    inline void copy(const CollSrc& src, CollDest& dest)
    {
	vampir_trace<3003> tracer;
	return copy(src, traits::flatcat2<CollSrc, tag::matrix, tag::vector>(),
		    dest, traits::flatcat2<CollDest, tag::matrix, tag::vector>());
    }


} // namespace mtl

#endif // MTL_COPY_INCLUDE
