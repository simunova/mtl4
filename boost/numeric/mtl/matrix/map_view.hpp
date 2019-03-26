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

#ifndef MTL_MAP_VIEW_INCLUDE
#define MTL_MAP_VIEW_INCLUDE

#include <utility>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/utility/is_multi_vector_expr.hpp>
#include <boost/numeric/mtl/matrix/crtp_base_matrix.hpp>
#include <boost/numeric/mtl/operation/sub_matrix.hpp>
#include <boost/numeric/mtl/operation/sfunctor.hpp>
#include <boost/numeric/mtl/operation/tfunctor.hpp>
#include <boost/numeric/mtl/operation/conj.hpp>
#include <boost/numeric/mtl/operation/imag.hpp>
#include <boost/numeric/mtl/operation/real.hpp>

#include <boost/numeric/mtl/matrix/mat_expr.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>


namespace mtl { namespace mat { namespace detail {
    // Forward declaration for friend declaration
    template <typename, typename> struct map_value;

    template <typename Functor, typename Matrix> struct map_vector {}; // not defined in general
    template <typename Functor, typename Vector>
    struct map_vector<Functor, mtl::mat::multi_vector<Vector> >
    {
	typedef mtl::vec::map_view<Functor, Vector> type;
    };

}}}

namespace mtl { namespace mat {

template <typename Functor, typename Matrix> 
struct map_view 
  : public const_crtp_base_matrix< map_view<Functor, Matrix>, 
				   typename Functor::result_type, typename Matrix::size_type >,
    public mat_expr< map_view<Functor, Matrix> >
{
    typedef map_view                                   self;
    typedef mat_expr< self >                           expr_base;
    typedef Matrix                                     other;
    typedef const Matrix&                              const_ref_type;
    typedef typename Matrix::orientation               orientation;
 
    typedef typename Functor::result_type              value_type;
    typedef typename Functor::result_type              const_reference;

    typedef typename Matrix::key_type                  key_type;
    typedef typename Matrix::size_type                 size_type;
    typedef typename Matrix::index_type                index_type;
    typedef typename Matrix::dim_type                  dim_type;

    struct dummy { typedef void type; };
    typedef typename boost::mpl::eval_if<mtl::traits::is_multi_vector_expr<Matrix>, detail::map_vector<Functor, Matrix>, dummy>::type vector_type;

    map_view (const Functor& functor, const other& ref) : functor(functor), ref(ref) {}
    
    map_view (const Functor& functor, boost::shared_ptr<Matrix> p) 
      : functor(functor), my_copy(p), ref(*p) {}
    
#ifdef MTL_WITH_MOVE    
  map_view (self&& that) : my_copy(std::move(that.my_copy)), functor(that.functor), ref(that.ref) {}
  map_view (const self& that) : functor(that.functor), ref(that.ref) { assert(that.my_copy.use_count() == 0); }
#endif


    value_type operator() (size_type r, size_type c) const
    { 
        return functor(ref(r, c));
    }
    // for multi_vector, needs enable_if since only defined for multi_vector
    template <typename S>
    typename boost::lazy_enable_if<boost::is_integral<S>, detail::map_vector<Functor, Matrix> >::type
    vector(S c) const
    {
	return typename detail::map_vector<Functor, Matrix>::type(functor, ref.vector(c));
    }

    size_type dim1() const { return ref.dim1(); }
    size_type dim2() const { return ref.dim2(); }
    dim_type dimensions() const { return ref.dimensions(); }

    size_type begin_row() const { return ref.begin_row(); }
    size_type end_row() const { return ref.end_row(); }
    size_type begin_col() const { return ref.begin_col(); }
    size_type end_col() const {	return ref.end_col(); }
    
    size_type nnz() const { return ref.nnz(); }

    friend size_type inline num_rows(const self& A) 
    { 	using mtl::mat::num_rows; return num_rows(A.ref);     }
    friend size_type inline num_cols(const self& A) 
    { 	using mtl::mat::num_cols; return num_cols(A.ref);     }
    template <typename, typename> friend struct detail::map_value;

  protected:
    boost::shared_ptr<Matrix>           my_copy;
  public:
    Functor           functor;
    const other&      ref;
};
   
template <typename Functor, typename Matrix> 
inline std::size_t size(const map_view<Functor, Matrix>& A)
{     return num_rows(A) * num_rows(A); }

// ==========
// Sub matrix
// ==========

template <typename Functor, typename Matrix>
struct sub_matrix_t< mtl::mat::map_view<Functor, Matrix> >
{
    typedef mtl::mat::map_view<Functor, Matrix>                                     view_type;

    // Mapping of sub-matrix type
    typedef typename sub_matrix_t<Matrix>::const_sub_matrix_type                       ref_sub_type;
    typedef mtl::mat::map_view<Functor, ref_sub_type>                               const_sub_matrix_type;
    typedef typename view_type::size_type                                              size_type;

    const_sub_matrix_type operator()(view_type const& view, size_type begin_r, size_type end_r, 
				     size_type begin_c, size_type end_c)
    {
	typedef boost::shared_ptr<ref_sub_type>                        pointer_type;

	// Submatrix of referred matrix (or view)
	// Create a submatrix, whos address will be kept by map_view
	// Functor is copied from view
	pointer_type p(new ref_sub_type(sub_matrix(view.ref, begin_r, end_r, begin_c, end_c)));
	return const_sub_matrix_type(view.functor, p); 
    }
};


}} // namespace mtl::matrix


namespace mtl { namespace traits {

    namespace detail {


	template <typename Functor, typename Matrix> 
	struct map_value
	{
	    typedef typename Matrix::key_type                      key_type;
	    typedef typename mtl::mat::map_view<Functor, Matrix>::value_type value_type;
    	
	    map_value(mtl::mat::map_view<Functor, Matrix> const& map_matrix) 
	      : map_matrix(map_matrix), its_value(map_matrix.ref) 
	    {}

	    value_type operator() (key_type const& key) const
	    {
		return map_matrix.functor(its_value(key));
	    }

	  protected:
	    mtl::mat::map_view<Functor, Matrix> const&        map_matrix;
		typename mtl::traits::const_value<Matrix>::type its_value;
        };



	template <typename Functor, typename Matrix> 
	struct mapped_row
	{
	    typedef typename Matrix::key_type   key_type;
	    typedef typename Matrix::size_type  size_type;
    	
		explicit mapped_row(const mtl::mat::map_view<Functor, Matrix>& view) : its_row(view.ref) {}
		explicit mapped_row(const mtl::mat::banded_view<Matrix>& view) : its_row(view.ref) {}

	    size_type operator() (key_type const& key) const
	    {
		return its_row(key);
	    }

	  protected:
	    typename row<Matrix>::type  its_row;
        };


        template <typename Functor, typename Matrix> 
        struct mapped_col
        {
	    typedef typename Matrix::key_type   key_type;
	    typedef typename Matrix::size_type  size_type;
    	
	    mapped_col(const mtl::mat::map_view<Functor, Matrix>& view) : its_col(view.ref) {}
	    mapped_col(const mtl::mat::banded_view<Matrix>& view) : its_col(view.ref) {}

	    size_type operator() (key_type const& key) const
	    {
		return its_col(key);
	    }

          protected:
	    typename col<Matrix>::type  its_col;
        };
	
    } // namespace detail
        
    template <typename Functor, typename Matrix> 
    struct row<mtl::mat::map_view<Functor, Matrix> >
    {
	typedef detail::mapped_row<Functor, Matrix>   type;
    };

    template <typename Functor, typename Matrix> 
    struct col<mtl::mat::map_view<Functor, Matrix> >
    {
	typedef detail::mapped_col<Functor, Matrix>   type;
    };

    template <typename Functor, typename Matrix> 
    struct const_value<mtl::mat::map_view<Functor, Matrix> >
    {
	typedef detail::map_value<Functor, Matrix>  type;
    };


    // ================
    // Range generators
    // ================

    // Use range_generator of original matrix
    template <typename Tag, typename Functor, typename Matrix> 
    struct range_generator<Tag, mtl::mat::map_view<Functor, Matrix> >
	: public detail::referred_range_generator<mtl::mat::map_view<Functor, Matrix>, 
						  range_generator<Tag, Matrix> >
    {};

    // To disambigue
    template <typename Functor, typename Matrix> 
    struct range_generator<tag::major, mtl::mat::map_view<Functor, Matrix> >
	: public detail::referred_range_generator<mtl::mat::map_view<Functor, Matrix>, 
						  range_generator<tag::major, Matrix> >
    {};


}} // mtl::traits


namespace mtl { namespace mat {

template <typename Scaling, typename Matrix>
struct scaled_view
  : public map_view<tfunctor::scale<Scaling, typename Matrix::value_type>, Matrix>
{
    typedef tfunctor::scale<Scaling, typename Matrix::value_type>  functor_type;
    typedef map_view<functor_type, Matrix>                         base;
    typedef scaled_view                                            self;

    scaled_view(const Scaling& scaling, const Matrix& matrix)
      : base(functor_type(scaling), matrix)
    {}
    
    scaled_view(const Scaling& scaling, boost::shared_ptr<Matrix> p)
      : base(functor_type(scaling), p)
    {}

#ifdef MTL_WITH_MOVE    
    scaled_view (self&& that) : base(that) {}
    scaled_view (const self& that) : base(that) {}
#endif
};

// rscaled_view -- added by Hui Li
template <typename Matrix, typename RScaling>
struct rscaled_view
  : public map_view<tfunctor::rscale<typename Matrix::value_type,RScaling>, Matrix>
{
    typedef tfunctor::rscale<typename Matrix::value_type, RScaling>  functor_type;
    typedef map_view<functor_type, Matrix>                          base;
    typedef rscaled_view                                            self;
	
    rscaled_view(const Matrix& matrix, const RScaling& rscaling)
      : base(functor_type(rscaling),matrix)
    {}

    rscaled_view(boost::shared_ptr<Matrix> p, const RScaling& rscaling)
      : base(functor_type(rscaling), p)
    {}

#ifdef MTL_WITH_MOVE    
    rscaled_view (self&& that) : base(that) {}
    rscaled_view (const self& that) : base(that) {}
#endif
};
	
// divide_by_view -- added by Hui Li
template <typename Matrix, typename Divisor>
struct divide_by_view
  : public map_view<tfunctor::divide_by<typename Matrix::value_type,Divisor>, Matrix>
{
    typedef tfunctor::divide_by<typename Matrix::value_type, Divisor>  functor_type;
    typedef map_view<functor_type, Matrix>                             base;
    typedef divide_by_view                                             self;
	
    divide_by_view(const Matrix& matrix,const Divisor& div)
      : base(functor_type(div), matrix)
    {}
	
    divide_by_view(boost::shared_ptr<Matrix> p, const Divisor& div)
      : base(functor_type(div), p)
    {}
	
#ifdef MTL_WITH_MOVE    
    divide_by_view (self&& that) : base(that) {}
    divide_by_view (const self& that) : base(that) {}
#endif
};

template <typename Matrix>
struct conj_view
  : public map_view<mtl::sfunctor::conj<typename Matrix::value_type>, Matrix>
{
    typedef mtl::sfunctor::conj<typename Matrix::value_type>            functor_type;
    typedef map_view<functor_type, Matrix>                              base;
    typedef conj_view                                                   self;

    conj_view(const Matrix& matrix) : base(functor_type(), matrix) {}
    conj_view(boost::shared_ptr<Matrix> p) : base(functor_type(), p) {}

#ifdef MTL_WITH_MOVE    
    conj_view (self&& that) : base(that) {}
    conj_view (const self& that) : base(that) {}
#endif
};

template <typename Matrix>
struct imag_view
  : public map_view<mtl::sfunctor::imag<typename Matrix::value_type>, Matrix>
{
    typedef mtl::sfunctor::imag<typename Matrix::value_type>            functor_type;
    typedef map_view<functor_type, Matrix>                              base;
    typedef imag_view                                                   self;

    imag_view(const Matrix& matrix) : base(functor_type(), matrix) {}
    imag_view(boost::shared_ptr<Matrix> p) : base(functor_type(), p) {}

#ifdef MTL_WITH_MOVE    
    imag_view (self&& that) : base(that) {}
    imag_view (const self& that) : base(that) {}
#endif
};

template <typename Matrix>
struct negate_view
  : public map_view<mtl::sfunctor::negate<typename Matrix::value_type>, Matrix>
{
    typedef mtl::sfunctor::negate<typename Matrix::value_type>            functor_type;
    typedef map_view<functor_type, Matrix>                              base;
    typedef negate_view                                                   self;

    negate_view(const Matrix& matrix) : base(functor_type(), matrix) {}
    negate_view(boost::shared_ptr<Matrix> p) : base(functor_type(), p) {}

#ifdef MTL_WITH_MOVE    
    negate_view (self&& that) : base(that) {}
    negate_view (const self& that) : base(that) {}
#endif
};

template <typename Matrix>
struct real_view
  : public map_view<mtl::sfunctor::real<typename Matrix::value_type>, Matrix>
{
    typedef mtl::sfunctor::real<typename Matrix::value_type>            functor_type;
    typedef map_view<functor_type, Matrix>                              base;
    typedef real_view                                                   self;

    real_view(const Matrix& matrix) : base(functor_type(), matrix) {}
    real_view(boost::shared_ptr<Matrix> p) : base(functor_type(), p) {}

#ifdef MTL_WITH_MOVE    
    real_view (self&& that) : base(that) {}
    real_view (const self& that) : base(that) {}
#endif
};

template <typename Matrix>
struct exp_view
  : public map_view<mtl::sfunctor::exp<typename Matrix::value_type>, Matrix>
{
    typedef mtl::sfunctor::exp<typename Matrix::value_type>            functor_type;
    typedef map_view<functor_type, Matrix>                              base;
    typedef exp_view                                                   self;

    exp_view(const Matrix& matrix) : base(functor_type(), matrix) {}
    exp_view(boost::shared_ptr<Matrix> p) : base(functor_type(), p) {}

#ifdef MTL_WITH_MOVE    
    exp_view (self&& that) : base(that) {}
    exp_view (const self& that) : base(that) {}
#endif
};





template <typename Scaling, typename Matrix>
struct sub_matrix_t< mtl::mat::scaled_view<Scaling, Matrix> >
  : public sub_matrix_t< mtl::mat::map_view<tfunctor::scale<Scaling, typename Matrix::value_type>, 
					       Matrix> >
{};

template <typename Matrix>
struct sub_matrix_t< mtl::mat::conj_view<Matrix> >
  : public sub_matrix_t< mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, Matrix> >
{};

template <typename Matrix, typename RScaling>
struct sub_matrix_t< mtl::mat::rscaled_view<Matrix, RScaling> >
  : public sub_matrix_t< mtl::mat::map_view<tfunctor::rscale<typename Matrix::value_type, RScaling>, 
					       Matrix> >
{};

template <typename Matrix, typename Divisor>
struct sub_matrix_t< mtl::mat::divide_by_view<Matrix, Divisor> >
  : public sub_matrix_t< mtl::mat::map_view<tfunctor::divide_by<typename Matrix::value_type, Divisor>, 
					       Matrix> >
{};


}} // namespace mtl::matrix

namespace mtl { namespace sfunctor {

    template <typename Matrix>
    struct conj_aux<Matrix, tag::matrix>
    {
	typedef mat::conj_view<Matrix> result_type;

	static inline result_type apply(const Matrix& matrix)
	{
	    return result_type(matrix);
	}

	result_type operator() (const Matrix& matrix) const
	{
	    return apply(matrix);
	}
    };

}} // namespace mtl::sfunctor

// Traits for specific views
namespace mtl { namespace traits {

template <typename Scaling, typename Matrix>
struct row< mtl::mat::scaled_view<Scaling, Matrix> >
  : public row< mtl::mat::map_view<tfunctor::scale<Scaling, typename Matrix::value_type>, 
				      Matrix> >
{};

template <typename Matrix>
struct row< mtl::mat::conj_view<Matrix> >
  : public row< mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, Matrix> >
{};

template <typename Matrix, typename RScaling>
struct row< mtl::mat::rscaled_view<Matrix, RScaling> >
  : public row< mtl::mat::map_view<tfunctor::rscale<typename Matrix::value_type, RScaling>, 
				      Matrix> >
{};

template <typename Matrix, typename Divisor>
struct row< mtl::mat::divide_by_view<Matrix, Divisor> >
  : public row< mtl::mat::map_view<tfunctor::divide_by<typename Matrix::value_type, Divisor>, 
				      Matrix> >
{};

template <typename Matrix>
struct row< mtl::mat::exp_view<Matrix> >
  : row< mtl::mat::map_view<mtl::sfunctor::exp<typename Matrix::value_type>, Matrix> >
{};


template <typename Scaling, typename Matrix>
struct col< mtl::mat::scaled_view<Scaling, Matrix> >
  : public col< mtl::mat::map_view<tfunctor::scale<Scaling, typename Matrix::value_type>, 
				      Matrix> >
{};

template <typename Matrix>
struct col< mtl::mat::conj_view<Matrix> >
  : public col< mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, Matrix> >
{};

template <typename Matrix, typename RScaling>
struct col< mtl::mat::rscaled_view<Matrix, RScaling> >
  : public col< mtl::mat::map_view<tfunctor::rscale<typename Matrix::value_type, RScaling>, 
				      Matrix> >
{};

template <typename Matrix, typename Divisor>
struct col< mtl::mat::divide_by_view<Matrix, Divisor> >
  : public col< mtl::mat::map_view<tfunctor::divide_by<typename Matrix::value_type, Divisor>, 
				      Matrix> >
{};

template <typename Matrix>
struct col< mtl::mat::exp_view<Matrix> >
  : col< mtl::mat::map_view<mtl::sfunctor::exp<typename Matrix::value_type>, Matrix> >
{};





template <typename Scaling, typename Matrix>
struct const_value< mtl::mat::scaled_view<Scaling, Matrix> >
  : public const_value< mtl::mat::map_view<tfunctor::scale<Scaling, typename Matrix::value_type>, 
					      Matrix> >
{};

template <typename Matrix>
struct const_value< mtl::mat::conj_view<Matrix> >
  : public const_value< mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, Matrix> >
{};

template <typename Matrix, typename RScaling>
struct const_value< mtl::mat::rscaled_view<Matrix, RScaling> >
  : public const_value< mtl::mat::map_view<tfunctor::rscale<typename Matrix::value_type, RScaling>, 
					      Matrix> >
{};

template <typename Matrix, typename Divisor>
struct const_value< mtl::mat::divide_by_view<Matrix, Divisor> >
  : public const_value< mtl::mat::map_view<tfunctor::divide_by<typename Matrix::value_type, Divisor>, 
					      Matrix> >
{};

template <typename Matrix>
struct const_value< mtl::mat::exp_view<Matrix> >
  : const_value< mtl::mat::map_view<mtl::sfunctor::exp<typename Matrix::value_type>, Matrix> >
{};




template <typename Tag, typename Scaling, typename Matrix>
struct range_generator< Tag, mtl::mat::scaled_view<Scaling, Matrix> >
    : public range_generator< Tag, mtl::mat::map_view<tfunctor::scale<Scaling, typename Matrix::value_type>, 
						    Matrix> >
{};

template <typename Tag, typename Matrix>
struct range_generator< Tag, mtl::mat::conj_view<Matrix> >
    : public range_generator< Tag, mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, Matrix> >
{};

template <typename Tag, typename Matrix, typename RScaling>
struct range_generator< Tag, mtl::mat::rscaled_view<Matrix, RScaling> >
    : public range_generator< Tag, mtl::mat::map_view<tfunctor::rscale<typename Matrix::value_type, RScaling>, 
						    Matrix> >
{};

template <typename Tag, typename Matrix, typename Divisor>
struct range_generator< Tag, mtl::mat::divide_by_view<Matrix, Divisor> >
    : public range_generator< Tag, mtl::mat::map_view<tfunctor::divide_by<typename Matrix::value_type, Divisor>, 
						    Matrix> >
{};

template <typename Tag, typename Matrix>
struct range_generator< Tag, mtl::mat::exp_view<Matrix> >
  : range_generator< Tag, mtl::mat::map_view<mtl::sfunctor::exp<typename Matrix::value_type>, Matrix> >
{};


template <typename Scaling, typename Matrix>
struct range_generator< tag::major, mtl::mat::scaled_view<Scaling, Matrix> >
    : public range_generator< tag::major, mtl::mat::map_view<tfunctor::scale<Scaling, typename Matrix::value_type>, 
						    Matrix> >
{};

template <typename Matrix>
struct range_generator< tag::major, mtl::mat::conj_view<Matrix> >
    : public range_generator< tag::major, mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, Matrix> >
{};

template <typename Matrix, typename RScaling>
struct range_generator< tag::major, mtl::mat::rscaled_view<Matrix, RScaling> >
    : public range_generator< tag::major, mtl::mat::map_view<tfunctor::rscale<typename Matrix::value_type, RScaling>, 
						    Matrix> >
{};

template <typename Matrix, typename Divisor>
struct range_generator< tag::major, mtl::mat::divide_by_view<Matrix, Divisor> >
    : public range_generator< tag::major, mtl::mat::map_view<tfunctor::divide_by<typename Matrix::value_type, Divisor>, 
						    Matrix> >
{};

template <typename Matrix>
struct range_generator< tag::major, mtl::mat::exp_view<Matrix> >
  : range_generator< tag::major, mtl::mat::map_view<mtl::sfunctor::exp<typename Matrix::value_type>, Matrix> >
{};


}} // mtl::traits


#endif // MTL_MAP_VIEW_INCLUDE
