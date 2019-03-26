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

#ifndef MTL_MATRIX_BANDED_VIEW_INCLUDE
#define MTL_MATRIX_BANDED_VIEW_INCLUDE

#include <utility>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/utility/parameters.hpp>
#include <boost/numeric/mtl/matrix/crtp_base_matrix.hpp>
#include <boost/numeric/mtl/matrix/base_matrix.hpp>
#include <boost/numeric/mtl/operation/sfunctor.hpp>
#include <boost/numeric/mtl/matrix/mat_expr.hpp>
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>


// Is not mutable because masking out some values forbids returning references
//
// Arbitrary combinations with other views (using shared_ptr) is planned

namespace mtl { namespace mat {

// Forward
namespace detail { 
    template <typename> struct banded_value; 
    template <typename, typename> struct mapped_value; 
}


template <typename Matrix> 
struct banded_view 
  : public const_crtp_base_matrix< banded_view<Matrix>, 
				   typename Matrix::value_type, typename Matrix::size_type >,
    public mat_expr< banded_view<Matrix> >,
    public base_matrix<typename Matrix::value_type, 
		       typename mtl::traits::parameters<Matrix>::type>
{
    typedef banded_view                                self;
    typedef mat_expr< self >                           expr_base;
    typedef typename mtl::traits::parameters<Matrix>::type parameters;

    typedef base_matrix<typename Matrix::value_type, parameters> base;
    
    typedef Matrix                                     other;
    typedef typename Matrix::orientation               orientation;
    typedef typename Matrix::index_type                index_type;
    // typedef typename Matrix::parameters                parameters;

    typedef typename Matrix::value_type                value_type;
    typedef typename Matrix::const_reference           const_reference;

    typedef typename Matrix::key_type                  key_type;
    typedef typename Matrix::size_type                 size_type;
    typedef typename Matrix::dim_type                  dim_type;

    typedef long int                                   bsize_type;

    banded_view(const other& ref, bsize_type begin, bsize_type end) 
      : base(dim_type(mtl::mat::num_rows(ref), mtl::mat::num_cols(ref)), ref.nnz()), 
	ref(ref), begin(begin), end(end) 
    {}

    banded_view(const boost::shared_ptr<Matrix>& p, bsize_type begin, bsize_type end) 
	: base(dim_type(mtl::mat::num_rows(*p), mtl::mat::num_cols(*p)), p->nnz()), 
	  my_copy(p), ref(*p), begin(begin), end(end) 
    {}

#ifdef MTL_WITH_MOVE    
    banded_view (self&& that) : my_copy(std::move(that.my_copy)), ref(that.ref), begin(that.begin), end(that.end) {}
    banded_view (const self& that) : ref(that.ref), begin(that.begin), end(that.end) { assert(that.my_copy.use_count() == 0); }
#endif

    value_type operator() (size_type r, size_type c) const
    {
	using math::zero;
	bsize_type bc= static_cast<bsize_type>(c), br= static_cast<bsize_type>(r),
	           band= bc - br;
	// Need value to return correct zero as well (i.e. matrices itself)
	value_type v= ref(r, c); 
	return begin <= band && band < end ? v : zero(v);
    }

    // need const functions
    bsize_type get_begin() const { return begin; }
    bsize_type get_end() const { return end; }

    template <typename> friend struct detail::banded_value;
    template <typename, typename> friend struct detail::map_value;
    //template <typename> friend struct ::mtl::sub_matrix_t<self>;

    friend size_type inline num_rows(const self& A) 
    { 	using mtl::mat::num_rows; return num_rows(A.ref);     }
    friend size_type inline num_cols(const self& A) 
    { 	using mtl::mat::num_cols; return num_cols(A.ref);     }

  protected:
    boost::shared_ptr<Matrix>           my_copy;
  public:
    const other&      ref;
    bsize_type        begin, end;
};

template <typename Matrix> 
inline std::size_t size(const banded_view<Matrix>& A)
{
    return num_rows(A) * num_rows(A); 
}

// ==========
// Sub matrix
// ==========

template <typename Matrix>
struct sub_matrix_t< mtl::mat::banded_view<Matrix> >
{
    typedef mtl::mat::banded_view<Matrix>                                           view_type;

    // Mapping of sub-matrix type
    typedef typename sub_matrix_t<Matrix>::sub_matrix_type                        ref_sub_type;
    typedef mtl::mat::banded_view<ref_sub_type>                                     const_sub_matrix_type;
    typedef mtl::mat::banded_view<ref_sub_type>                                     sub_matrix_type;
    typedef typename view_type::size_type                                         size_type;

    sub_matrix_type operator()(view_type const& view, size_type begin_r, size_type end_r, 
				     size_type begin_c, size_type end_c)
    {
	typedef boost::shared_ptr<ref_sub_type>                        pointer_type;

	// Submatrix of referred matrix (or view)
	// Create a submatrix, whos address will be kept by banded_view
	pointer_type p(new ref_sub_type(sub_matrix(view.ref, begin_r, end_r, begin_c, end_c)));
	return sub_matrix_type(p, view.begin, view.end); 
    }
};


}} // namespace mtl::matrix




namespace mtl { namespace traits {

    using mtl::mat::banded_view;

    template <typename Matrix> 
    struct row<banded_view<Matrix> >
    {
	// from map_view
	typedef detail::mapped_row<sfunctor::identity<typename Matrix::value_type>, Matrix>   type;
    };

    template <typename Matrix> 
    struct col<banded_view<Matrix> >
    {
	// from map_view
	typedef detail::mapped_col<sfunctor::identity<typename Matrix::value_type>, Matrix>   type;
    };

    namespace detail {

	template <typename Matrix> 
	struct banded_value
	{
	    typedef typename Matrix::key_type              key_type;
	    typedef typename Matrix::value_type            value_type;
	    typedef banded_view<Matrix>                    view_type;
    	
	    banded_value(view_type const& view) 
		: view(view), its_row(view.ref), its_col(view.ref), its_value(view.ref) 
	    {}

	    value_type operator() (key_type const& key) const
	    {
		using math::zero;
		typedef typename view_type::bsize_type   bsize_type;

		bsize_type br= static_cast<bsize_type>(its_row(key)), 
                           bc= static_cast<bsize_type>(its_col(key)),
		           band= bc - br;
		// Need value to return correct zero as well (i.e. matrices itself)
		const value_type v= its_value(key);

		return view.get_begin() <= band && band < view.get_end() ? v : zero(v);
	    }

	  protected:
	    view_type const&                                view;
	    typename row<Matrix>::type         its_row;
	    typename col<Matrix>::type         its_col;
	    typename const_value<Matrix>::type its_value;
        };

    } // detail

    template <typename Matrix> 
    struct const_value<banded_view<Matrix> >
    {
	typedef detail::banded_value<Matrix>  type;
    };

    // ================
    // Range generators
    // ================

    // Use range_generator of original matrix
    template <typename Tag, typename Matrix> 
    struct range_generator<Tag, banded_view<Matrix> >
	: public detail::referred_range_generator<banded_view<Matrix>, range_generator<Tag, Matrix> >
    {};

#if 0 // It is more complicated than this because referred_range_generator returns Matrix's cursor and we 
      // cannot dispatch on this anymore
    template <typename Matrix> 
    struct range_generator<glas::tag::nz, 
			   typename detail::referred_range_generator<banded_view<Matrix>, range_generator<Tag, Matrix> >::type> 

			   detail::sub_matrix_cursor<banded_view<Matrix>, glas::tag::row, 2> >
    {
	typedef range_generator<glas::tag::row, banded_view<Matrix> >                    collection_type;
	typedef range_generator<glas::tag::nz, range_generator<glas::tag::row, Matrix> > other_generator;
	static int const                            level = other_generator::level;
	typedef typename other_generator::type       type;

	type begin(const collection_type& c)
	{
	    return 

	// 
	typedef typename range_generator<glas::tag::nz, detail::sub_matrix_cursor<Matrix>, glas::tag::row, 2>::type type;
#endif


    // To disambiguate
    template <typename Matrix> 
    struct range_generator<tag::major, banded_view<Matrix> >
	: public detail::referred_range_generator<banded_view<Matrix>, range_generator<tag::major, Matrix> >
    {};






}} // mtl::traits


#endif // MTL_MATRIX_BANDED_VIEW_INCLUDE
